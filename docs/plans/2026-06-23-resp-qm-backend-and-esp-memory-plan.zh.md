# XpongeCPP 与 Xponge-origin RESP/QM/ESP 内存治理实施方案

## 1. 目标

本方案覆盖两个直接目标：

1. 让 `Xponge-origin` 的 `RESP` 走统一 `QM` 后端，而不是继续绑死 `PySCF`
2. 让 `XpongeCPP` 与 `Xponge-origin` 的 `RESP ESP` 计算都支持受内存预算控制的自动分块，避免 `aux_e2` 全量中间张量导致 OOM

本方案同时要求：

- `RESP` 的用户入口保持兼容
- `PySCF` 与 `Psi4` 的后端选择语义保持统一
- `MCPB` 的本地 `RESP` 重拟合能自动受益

## 2. 非目标

本轮不包含：

- 重写 `SCF`、积分、`ESP` 或 `Hessian` 的数值内核
- 把 `RESP` 全部迁进 `C++`
- 一次性重构所有量化学功能
- 对所有后端强行实现同一套底层分块算法

## 3. 当前现状

### 3.1 XpongeCPP

- `RESP` 编排层已经存在统一后端入口：
  - [XpongeCPP/assign/resp.py](/mnt/data8t/Software/Xponge/Xponge-CPP/src/XpongeCPP/assign/resp.py)
- `QM` 抽象层已经存在：
  - [XpongeCPP/qm/models.py](/mnt/data8t/Software/Xponge/Xponge-CPP/src/XpongeCPP/qm/models.py)
  - [XpongeCPP/qm/scheduler.py](/mnt/data8t/Software/Xponge/Xponge-CPP/src/XpongeCPP/qm/scheduler.py)
  - [XpongeCPP/qm/backends/pyscf_backend.py](/mnt/data8t/Software/Xponge/Xponge-CPP/src/XpongeCPP/qm/backends/pyscf_backend.py)
  - [XpongeCPP/qm/backends/psi4_backend.py](/mnt/data8t/Software/Xponge/Xponge-CPP/src/XpongeCPP/qm/backends/psi4_backend.py)
- `RESP` 数值核心已经与 `QM` 后端解耦，`C++ core` 只负责 `MK grid` 和 `RESP fit`：
  - [XpongeCPP/assign/resp_core.py](/mnt/data8t/Software/Xponge/Xponge-CPP/src/XpongeCPP/assign/resp_core.py)

已确认问题：

- `PySCF` 后端当前仍默认走全量 `aux_e2`：
  - [compute_esp()](/mnt/data8t/Software/Xponge/Xponge-CPP/src/XpongeCPP/qm/backends/pyscf_backend.py:107)
- `Psi4` 后端当前也默认把全部 grid 一次性交给 `compute_esp_over_grid_in_memory(...)`：
  - [compute_esp()](/mnt/data8t/Software/Xponge/Xponge-CPP/src/XpongeCPP/qm/backends/psi4_backend.py:90)

结论：

- `XpongeCPP` 的 `RESP` 已经走统一 `QM` 后端
- `ESP` 内存问题仍然存在
- 最合适的修复点是 `qm/backends/*_backend.py`，不是 `RESP core`

### 3.2 Xponge-origin

- `RESP` 仍是单体式 `PySCF` 实现：
  - [Xponge/assign/resp.py](/mnt/data8t/Software/Xponge/Xponge-origin/Xponge/assign/resp.py)
- `Assign.calculate_charge("RESP")` 直接调用它：
  - [Xponge/assign/__init__.py](/mnt/data8t/Software/Xponge/Xponge-origin/Xponge/assign/__init__.py:751)

已确认问题：

- 顶层直接 `import pyscf`
- 没有统一 `QM` 后端抽象
- `ESP` 计算沿用全量 `aux_e2` 再 `MemoryError` fallback 的策略
- 无法让 `RESP`、未来 `MCPB`、未来 `Hessian` 共用同一后端层

### 3.3 已知 OOM 证据

对 `G6141.mol2` 的实际检查结果为：

- `natm = 53`
- `NAO = 538`
- 默认 `RESP` grid 数约 `35382`
- `df.incore.aux_e2(mol, fakemol)` 尝试分配形状近似为：
  - `(538, 538, 35382, 1)`
- 对应内存约：
  - `76.3 GiB`

因此当前的主要 OOM 点不是：

- `SCF`
- `matrix_a0`
- `RESP` 线性求解

而是：

- `ESP/MEP` 生成阶段的一次性全量中间张量

## 4. 设计结论

### 4.1 统一接口结论

两个仓库都应统一采用：

- `RESP` 只负责编排
- `QM` 后端负责：
  - `SCF`
  - `geometry optimization`
  - `ESP on grid`
  - `Hessian`

`RESP` 不再直接了解 `PySCF` / `Psi4` 的分子构建和积分细节。

### 4.2 内存治理结论

用户应暴露的是：

- `ESP` 内存预算

而不是：

- 手工指定 `AO` 分块数
- 手工指定 `grid` 分块数

建议的用户参数：

```python
assign.calculate_charge(
    "resp",
    backend="pyscf",
    esp_memory_limit="1GB",
)
```

内部再由后端自动选择：

- `full`
- `grid_chunk`
- `shell_grid_chunk`
- `pointwise`

### 4.3 默认策略结论

若用户未显式给出预算，建议默认：

- `esp_memory_limit = 1 GiB`

理由：

- 对 `RESP` 这种辅助型 `ESP` 计算来说，`1 GiB` 已足够大
- 可显著降低“默认调用直接炸机器”的风险
- 仍保留高级用户手工放宽预算的能力

## 5. 统一接口设计

## 5.1 用户层接口

两个仓库最终统一支持：

```python
assign.calculate_charge(
    "resp",
    backend="pyscf",
    basis="6-31g*",
    esp_memory_limit="1GB",
    esp_chunk_policy="auto",
)
```

参数含义：

- `backend`
  - `pyscf` / `psi4`
- `esp_memory_limit`
  - 预算上限，支持字符串或整数 bytes
- `esp_chunk_policy`
  - `auto`
  - `full`
  - `grid`
  - `dual`
  - `pointwise`

其中：

- `auto` 为默认策略
- `dual` 表示 `AO(shell) + grid` 双向分块

## 5.2 QM 请求模型

建议把 `ESP` 专属配置放到 `ESPGridRequest`，而不是强塞进 `QMRunOptions.memory`。

目标扩展字段：

```python
ESPGridRequest(
    grid_points_bohr=...,
    include_nuclear_term=False,
    include_electronic_term=True,
    memory_limit_bytes=1073741824,
    chunk_policy="auto",
    safety_factor=0.8,
)
```

原因：

- `SCF memory` 与 `ESP chunk memory` 是两个不同层面的约束
- `RESP` 关心的是 `ESP` 阶段峰值，不是一般 `SCF` 内存

## 5.3 ESP 结果诊断

建议在 `ESPResult` 中新增诊断元数据，例如：

- `mode`
  - `full`
  - `grid_chunk`
  - `shell_grid_chunk`
  - `pointwise`
- `estimated_full_bytes`
- `memory_limit_bytes`
- `grid_chunk_count`
- `shell_block_count`

这样：

- 调试方便
- 测试可验证
- 用户日志可解释

## 6. 分块算法设计

## 6.1 共同原则

最终需要的结果只是：

- 每个采样点的 `electronic_esp[p]`

因此不应默认长期持有：

\[
T_{ijp} \sim AO \times AO \times N_{grid}
\]

的完整中间张量。

应采用：

- 先估算内存
- 在预算内优先选吞吐更高的路径
- 超预算时自动切块

## 6.2 PySCF 后端策略

### 路径 A：full

当估算的全量张量大小满足：

\[
\text{estimated\_bytes} \le \text{memory\_limit} \times \text{safety\_factor}
\]

时，保留现有：

- `fakemol_for_charges(...)`
- `df.incore.aux_e2(...)`
- `einsum(...)`

### 路径 B：grid chunk

当全量超预算，但单个完整 `AO × AO × Nchunk` 仍能放进预算时：

- 把 `grid_points` 分块
- 每次对一个 grid chunk 构造 `fakemol`
- 对该块调用 `aux_e2`
- 立刻收缩成该块 `electronic_esp`
- 写入输出向量后释放中间块

该路径是第一阶段应优先交付的主修复。

### 路径 C：shell + grid dual chunk

当 `grid chunk` 仍然不够稳，或未来 `NAO` 更大时，增加 `PySCF shell slice` 双向分块。

已确认 `PySCF` 接口支持：

```python
df.incore.aux_e2(..., shls_slice=(ish0, ish1, jsh0, jsh1, aux0, aux1))
```

因此可采用：

- 外层切 `aux/grid` 块
- 内层切 `AO shell-pair` 块
- 用 `mol.ao_loc_nr()` 把 shell block 映射到 AO block
- 每一块立即：
  - 取出 `dm[iao0:iao1, jao0:jao1]`
  - 计算局部 `aux_e2`
  - 累加到当前 chunk 的 `electronic_esp`

注意：

- `dual chunk` 实际上应以 `shell` 为切分单元，而不是裸 `AO index`
- 第一版可固定使用 `aosym='s1'`，先保证正确性

### 路径 D：pointwise fallback

若预算极小，或 `shell/grid` 分块后仍难以稳定满足预算，则回退到：

- 单点 `int1e_rinv`
- 与密度矩阵收缩

这是最慢但最稳的兜底路径。

## 6.3 Psi4 后端策略

`Psi4` 当前接口没有像 `PySCF aux_e2` 那样暴露 `AO/shell` 级切块入口。  
因此本轮对 `Psi4` 的策略应为：

- 支持 `grid chunk`
- 不承诺第一版实现 `AO` 方向分块

即：

- 把 `grid_points_bohr` 分成多块
- 反复调用 `compute_esp_over_grid_in_memory(...)`
- 拼接结果

这样虽然不具备 `dual chunk`，但接口仍然统一，后端内部策略允许不同。

## 6.4 分阶段实现建议

为降低风险，建议分两阶段交付：

### Stage A：先解决当前 OOM

- `auto`
- `full`
- `grid_chunk`
- `pointwise`

### Stage B：再上 shell/grid dual chunk

- 仅对 `PySCF` 后端实现
- 主要面向更大 `NAO` 体系

这样可以先快速修掉当前 `G6141` 级别问题，再补更大体系的可扩展性。

## 7. XpongeCPP 实施计划

## 7.1 Phase 1：扩展 `QM` 请求与结果模型

修改文件：

- [qm/models.py](/mnt/data8t/Software/Xponge/Xponge-CPP/src/XpongeCPP/qm/models.py)
- [qm/scheduler.py](/mnt/data8t/Software/Xponge/Xponge-CPP/src/XpongeCPP/qm/scheduler.py)

计划改动：

- 为 `ESPGridRequest` 增加：
  - `memory_limit_bytes`
  - `chunk_policy`
  - `safety_factor`
- 为 `ESPResult` 增加：
  - `diagnostics`
- 让 `compute_esp_on_grid(...)` 支持透传这些配置

## 7.2 Phase 2：改造 `PySCF` `compute_esp()`

修改文件：

- [qm/backends/pyscf_backend.py](/mnt/data8t/Software/Xponge/Xponge-CPP/src/XpongeCPP/qm/backends/pyscf_backend.py)

计划改动：

- 抽出 helper：
  - `_estimate_aux_e2_bytes(...)`
  - `_compute_esp_full(...)`
  - `_compute_esp_grid_chunked(...)`
  - `_compute_esp_shell_grid_chunked(...)`
  - `_compute_esp_pointwise(...)`
- `auto` 策略按预算自动选择路径
- 在 `ESPResult.diagnostics` 中记录实际选择

说明：

- `RESP core` 不需要改
- `MCPB local RESP refit` 会自动受益

## 7.3 Phase 3：改造 `Psi4` `compute_esp()`

修改文件：

- [qm/backends/psi4_backend.py](/mnt/data8t/Software/Xponge/Xponge-CPP/src/XpongeCPP/qm/backends/psi4_backend.py)

计划改动：

- 加入 `grid chunk` 支持
- 透传并记录相同的预算/策略字段
- 明确 `dual chunk` 当前仅对 `PySCF` 开放

## 7.4 Phase 4：把新参数接到 `RESP` 用户入口

修改文件：

- [assign/resp.py](/mnt/data8t/Software/Xponge/Xponge-CPP/src/XpongeCPP/assign/resp.py)
- 必要时同步兼容层：
  - [XpongeCPP/_compat/assign.py](/mnt/data8t/Software/Xponge/Xponge-CPP/src/XpongeCPP/_compat/assign.py)

计划改动：

- `resp_fit(...)` 接受：
  - `esp_memory_limit`
  - `esp_chunk_policy`
- 转换为 `ESPGridRequest` 的字段

## 7.5 Phase 5：测试与文档

测试重点：

- estimator 单测
- `tiny budget` 强制走 `grid chunk`
- `PySCF` shell/grid dual chunk 路径的单测
- `Psi4` 轻量 `grid chunk` smoke test
- 现有 `RESP` / `MCPB` 轻量回归

建议新增：

- repo 内小分子真实用例
- repo 外 `G6141.mol2` 脚本型 benchmark，不进 CI

文档更新：

- `RESP` 使用说明
- `QM backend` 使用说明
- `Windows/Psi4` 提示不变

## 8. Xponge-origin 实施计划

## 8.1 分支原则

`Xponge-origin` 的实现应在：

- `metal_assignment`

分支上继续推进，而不是直接落在 `master`。

## 8.2 Phase 1：迁回统一 `QM` 架构

目标：

- 让 `origin` 拥有与 `XpongeCPP` 对齐的：
  - `qm/models.py`
  - `qm/scheduler.py`
  - `qm/backends/pyscf_backend.py`
  - `qm/backends/psi4_backend.py`

建议文件：

- 新增 [Xponge/qm/](/mnt/data8t/Software/Xponge/Xponge-origin/Xponge/qm)
- 参考 `XpongeCPP` 当前实现迁移

## 8.3 Phase 2：重构 `RESP`

目标：

- 保留 `Assign.calculate_charge("RESP", ...)` 入口不变
- 让 `RESP` 改成编排层 + 数值核心两层

建议结构：

- [Xponge/assign/resp.py](/mnt/data8t/Software/Xponge/Xponge-origin/Xponge/assign/resp.py)
  - 编排层
- 新增 `Xponge/assign/resp_core.py`
  - 保留现有 Python 数值逻辑

说明：

- `origin` 不必把 `RESP fit core` 也立即迁成 `C++`
- 本轮重点是 `QM backend` 解耦与 `ESP` 内存治理

## 8.4 Phase 3：移植 `ESP` 预算接口

目标：

- `origin RESP` 公开支持：
  - `backend`
  - `esp_memory_limit`
  - `esp_chunk_policy`

实现方式：

- 优先复用 `XpongeCPP` 已验证过的后端逻辑
- 不在老 `master` 单体式 `resp.py` 上重复造一遍

## 8.5 Phase 4：测试

测试重点：

- `RESP` 默认 `PySCF` 路径不回归
- `Psi4` 轻量 capability / smoke test
- `tiny budget` 下触发 chunk 路径
- `G6141` 风格大 grid 场景的脚本验证

## 9. 推荐实施顺序

推荐顺序必须是：

1. 先修 `XpongeCPP`
2. 再把 `QM + ESP chunking` 迁回 `Xponge-origin`

原因：

- `XpongeCPP` 现在已经有统一 `QM` 层，修复点最集中
- `origin` 还没有这层结构
- 若先在 `origin` 老 `resp.py` 上单独修，会做两遍类似工作

## 10. 验收标准

两个仓库都应满足：

1. `RESP` 用户入口支持统一后端参数
2. `PySCF` 默认不再先赌全量 `aux_e2` 再靠 `MemoryError` 救场
3. 在 `esp_memory_limit="1GB"` 下，`G6141` 级别体系不再尝试分配 `76 GiB` 级中间张量
4. `MCPB local RESP refit` 自动复用改进后的 `ESP` 计算路径
5. `Psi4` 路径至少支持统一接口下的 `grid chunk`

## 11. 建议的提交切分

### XpongeCPP

1. `Add ESP memory-budget request model`
2. `Chunk PySCF RESP ESP evaluation`
3. `Chunk Psi4 RESP ESP evaluation`
4. `Expose RESP esp_memory_limit options`
5. `Add RESP/QM chunking tests and docs`

### Xponge-origin

1. `Port unified QM backend skeleton`
2. `Split RESP orchestration from numerical core`
3. `Port ESP memory-budgeted backend evaluation`
4. `Expose RESP backend and esp_memory_limit options`
5. `Add origin RESP/QM chunking tests`

## 12. 最终建议

这件事的最佳落点不是 `RESP core`，而是：

- `XpongeCPP`：`qm/backends/*_backend.py`
- `Xponge-origin`：先迁 `qm/`，再在相同层面落地

用户暴露 `1 GiB` 这种预算式接口是合理的。  
但工程实现上应遵循：

- 第一阶段先交付 `grid chunk`
- 第二阶段再交付 `PySCF shell/grid dual chunk`

这样能先可靠解决当前 OOM，再把设计扩展到更大 `NAO` 体系。
