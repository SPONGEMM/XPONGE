# Xponge Amber GLYCAM 与糖蛋白模板对齐修订计划

## 1. 目标

本计划的目标是让 `Xponge-origin` 在 Amber/GLYCAM 糖与糖蛋白模板支持上，先补齐当前最影响可用性的缺口，再逐步逼近 `AmberTools26` 的模板覆盖边界。

本轮修订的直接目标分为三层：

1. 修复 `ff14SB` / `ff19SB` 中 `HYP` 相关的标准蛋白残基派生与端基映射缺口
2. 使 `Xponge` 的核心糖蛋白残基组织与 `AmberTools26 leaprc.GLYCAM_06j-1` 的工作流语义对齐
3. 为后续补齐 `GLYCAM_06j` 扩展糖模板建立清晰的覆盖审计、优先级和迁移路径

## 2. 范围

### 2.1 本次包含

- 检查并修正 `ff14SB` / `ff19SB` 的 `HYP`、`CHYP`、`NHYP` 支持边界
- 对齐 `glycoprotein` 相关核心模板：
  - `NLN`
  - `OLS`
  - `OLT`
  - `OLP`
  - `HYP`
  - `NHYP`
  - `CHYP`
- 建立 `AmberTools26 -> Xponge` 的 `GLYCAM` 模板覆盖审计机制
- 规划 `SO3`、`MEX`、`TBT`、`CA2`、`0TV`、`0AE`、`0AF`、`0GL` 等扩展模板的补齐顺序
- 补齐必要的回归测试和说明文档

### 2.2 本次不包含

- 一次性手工补全 `AmberTools26 GLYCAM_06j-1.prep` 的全部缺项
- 重写 `load_pdb()` 的糖化学自动识别逻辑
- 自动把普通 `ASN/SER/THR/HYP` 改写成 `NLN/OLS/OLT/OLP` 的通用推断器
- 立即实现与 `AmberTools26` 百分之百的名字级完全镜像

## 3. 当前现状

## 3.1 Xponge 当前状态

### 标准 Amber 蛋白模板

Xponge 当前的 `Amber ff14SB / ff19SB` 主体模板位于：

- [Xponge/forcefield/amber/ff14SB.mol2](/mnt/data8t/Software/Xponge/Xponge-origin/Xponge/forcefield/amber/ff14SB.mol2)
- [Xponge/forcefield/amber/ff19SB.mol2](/mnt/data8t/Software/Xponge/Xponge-origin/Xponge/forcefield/amber/ff19SB.mol2)

加载与头尾派生逻辑位于：

- [Xponge/forcefield/amber/ff14sb.py](/mnt/data8t/Software/Xponge/Xponge-origin/Xponge/forcefield/amber/ff14sb.py)
- [Xponge/forcefield/amber/ff19sb.py](/mnt/data8t/Software/Xponge/Xponge-origin/Xponge/forcefield/amber/ff19sb.py)

已确认：

- `ff14SB.mol2` 中存在 `HYP` 和 `CHYP`
- `ff19SB.mol2` 中存在 `HYP` 和 `CHYP`
- 两个 Python loader 的标准蛋白 `residues = ...` 派生列表都没有 `HYP`
- `ff19SB` 当前未见 `NHYP` 模板

### GLYCAM 与糖蛋白模板

Xponge 当前 GLYCAM 入口位于：

- [Xponge/forcefield/amber/glycam_06j/__init__.py](/mnt/data8t/Software/Xponge/Xponge-origin/Xponge/forcefield/amber/glycam_06j/__init__.py)

糖模板组织位于：

- [d_pyranose.py](/mnt/data8t/Software/Xponge/Xponge-origin/Xponge/forcefield/amber/glycam_06j/d_pyranose.py)
- [d_furanose.py](/mnt/data8t/Software/Xponge/Xponge-origin/Xponge/forcefield/amber/glycam_06j/d_furanose.py)
- `l_pyranose.py`
- `l_furanose.py`

糖蛋白桥接模板位于：

- [glycoprotein.py](/mnt/data8t/Software/Xponge/Xponge-origin/Xponge/forcefield/amber/glycam_06j/glycoprotein.py)

当前已确认：

- `glycoprotein.py` 已显式接入 `OLP`、`OLT`、`NLN`、`OLS`
- 并派生其 `N...` / `C...` 形式
- `terminal.mol2` 当前只含 `OME` 与 `ROH`
- 尚未形成与 Amber `GLYCAM_amino*.lib + leaprc.GLYCAM_06j-1` 同等完整的糖蛋白模板组织

## 3.2 AmberTools26 当前状态

AmberTools26 的相关入口与模板主要位于：

- [leaprc.protein.ff14SB](/mnt/data8t/Software/AmberTools26/ambertools26_src/dat/leap/cmd/leaprc.protein.ff14SB)
- [leaprc.protein.ff19SB](/mnt/data8t/Software/AmberTools26/ambertools26_src/dat/leap/cmd/leaprc.protein.ff19SB)
- [leaprc.GLYCAM_06j-1](/mnt/data8t/Software/AmberTools26/ambertools26_src/dat/leap/cmd/leaprc.GLYCAM_06j-1)
- [GLYCAM_06j-1.prep](/mnt/data8t/Software/AmberTools26/ambertools26_src/dat/leap/prep/GLYCAM_06j-1.prep)
- [GLYCAM_amino_06j_12SB.lib](/mnt/data8t/Software/AmberTools26/ambertools26_src/dat/leap/lib/GLYCAM_amino_06j_12SB.lib)
- [GLYCAM_aminont_06j_12SB.lib](/mnt/data8t/Software/AmberTools26/ambertools26_src/dat/leap/lib/GLYCAM_aminont_06j_12SB.lib)
- [GLYCAM_aminoct_06j_12SB.lib](/mnt/data8t/Software/AmberTools26/ambertools26_src/dat/leap/lib/GLYCAM_aminoct_06j_12SB.lib)

已确认：

- `ff14SB` 使用 `HYP -> HYP / CHYP`
- `ff19SB` 使用 `HYP -> NHYP / CHYP`
- `GLYCAM_amino*.lib` 提供：
  - `HYP`, `NLN`, `OLP`, `OLS`, `OLT`
  - `NHYP`, `NNLN`, `NOLP`, `NOLS`, `NOLT`
  - `CHYP`, `CNLN`, `COLP`, `COLS`, `COLT`
- `GLYCAM_06j-1.prep` 不仅包含基础糖单元，也包含：
  - 糖修饰基团
  - 糖封端/帽基
  - 特殊糖模板
  - 少量辅助单元和离子

## 3.3 关键差异摘要

当前差异可以归纳为两类：

1. 标准蛋白层的 `HYP` 支持没有完全与 Amber 的 `ff14SB / ff19SB` 对齐
2. `GLYCAM` 扩展模板层没有完整迁入 Xponge

其中：

- `ff14SB` 的问题主要是 `HYP` 没进入派生/映射逻辑
- `ff19SB` 的问题是 `HYP` 没进入派生/映射逻辑，且缺 `NHYP`
- `GLYCAM` 侧则缺一批功能团与特殊糖模板，如：
  - `SO3`
  - `MEX`
  - `TBT`
  - `CA2`
  - `0TV`
  - `0AE`
  - `0AF`
  - `0GL`
  - 以及同类扩展族

## 4. 设计原则

### 4.1 先补核心工作流，再扩展覆盖率

先解决用户最容易遇到的两个问题：

- `HYP` 在标准蛋白和糖蛋白语境下的端基支持
- 常见 glycoprotein 核心模板的稳定加载与连接

再推进大批量 `GLYCAM` 扩展模板覆盖。

### 4.2 不用单一通用循环硬套 `HYP`

`HYP` 不是普通残基，头部约束与 `PRO` 类似。  
因此不应仅仅把 `HYP` 塞进现有 `residues = ...` 列表，尤其是 `ff14SB` 场景下，因为：

- `ff14SB` 不要求 `NHYP`
- `HYP` 的头部几何条件也不应走普通残基分支

### 4.3 Phase 3 不是“只补糖封端”

扩展模板层需要区分三类对象：

1. 小功能团 / 封端块  
   例如 `SO3`、`MEX`、`OME`、`ROH`、`TBT`
2. 完整修饰糖模板  
   例如 `0TV`、`0AE`、`0AF`、`0GL`
3. 辅助单元  
   例如 `CA2`

因此 Phase 3 应以“扩展模板层审计与补齐”为定义，而不是简单叫“补封端”。

### 4.4 先建立审计脚本，再扩模板

如果不先建立 `Amber -> Xponge` 的覆盖审计脚本，后续补模板会变成零散手工维护，难以持续回归。

## 5. 三阶段实施方案

## 5.1 Phase 1：标准蛋白 `HYP` 对齐

### 目标

让 `Xponge` 的 `ff14SB / ff19SB` 在 `HYP`、`CHYP`、`NHYP` 的端基与映射语义上对齐 Amber。

### 文件范围

- [Xponge/forcefield/amber/ff14sb.py](/mnt/data8t/Software/Xponge/Xponge-origin/Xponge/forcefield/amber/ff14sb.py)
- [Xponge/forcefield/amber/ff19sb.py](/mnt/data8t/Software/Xponge/Xponge-origin/Xponge/forcefield/amber/ff19sb.py)
- [Xponge/forcefield/amber/ff19SB.mol2](/mnt/data8t/Software/Xponge/Xponge-origin/Xponge/forcefield/amber/ff19SB.mol2)

### 改动方案

#### Phase 1A：ff14SB

- 保持 `HYP` / `CHYP` 使用现有 `mol2` 模板
- 不把 `HYP` 简单加入现有通用 `residues` 循环
- 为 `HYP` 单独加显式头尾逻辑：
  - 基础残基：`HYP`
  - C 端：`CHYP`
  - 不要求 `NHYP`
- 头部 link 条件按 `PRO` 风格处理
- 增加 `PDB` 残基映射：
  - `HYP -> HYP / CHYP`

#### Phase 1B：ff19SB

- 补入 `NHYP` 模板
- 为 `HYP/NHYP/CHYP` 单独加显式头尾逻辑
- 增加 `PDB` 残基映射：
  - `HYP -> NHYP / CHYP`
- 同样按 `PRO` 风格处理头部条件

### 验收标准

- `ff14SB`：
  - 中间 `HYP` 可正常加载
  - C 端 `HYP` 正确映射为 `CHYP`
- `ff19SB`：
  - 中间 `HYP` 可正常加载
  - N 端 `HYP` 正确映射为 `NHYP`
  - C 端 `HYP` 正确映射为 `CHYP`

## 5.2 Phase 2：糖蛋白核心模板与工作流对齐

### 目标

让 `Xponge` 当前常见糖蛋白工作流与 `AmberTools26 leaprc.GLYCAM_06j-1` 的核心语义对齐。

### 文件范围

- [Xponge/forcefield/amber/glycam_06j/glycoprotein.py](/mnt/data8t/Software/Xponge/Xponge-origin/Xponge/forcefield/amber/glycam_06j/glycoprotein.py)
- `Xponge/forcefield/amber/glycam_06j/*.mol2`
- 可能新增测试文件

### 改动方案

#### Phase 2A：明确核心残基集

将糖蛋白核心支持定义为：

- `NLN`
- `OLS`
- `OLT`
- `OLP`
- `HYP`
- `NHYP`
- `CHYP`

其中：

- `NLN/OLS/OLT/OLP` 及其 `N.../C...` 形式由 `glycoprotein` 层负责
- `HYP/NHYP/CHYP` 由标准蛋白 FF 层负责
- 两者通过统一加载与测试协同，而不是重复定义

#### Phase 2B：收口加载约定

为用户建立清晰约定：

- `ff14sb + glycam_06j + glycoprotein`
- `ff19sb + glycam_06j + glycoprotein`

并在文档或模块说明中明确：

- 哪些残基来自蛋白 FF
- 哪些残基来自 GLYCAM glycoprotein
- 依赖 `LINK` 的地方有哪些

#### Phase 2C：补核心回归测试

至少覆盖：

- `NLN`
- `OLS`
- `OLT`
- `OLP`
- `HYP/CHYP`
- `HYP/NHYP/CHYP`（ff19SB）

测试目标：

- 模板可解析
- 头尾映射正确
- 蛋白-糖连接位点可建立

### 验收标准

- `ff14SB + GLYCAM` 下，常见 `NLN/OLS/OLT/OLP` 糖蛋白路径可工作
- `ff19SB + GLYCAM` 下，`NHYP` 路径可工作
- 用户能明确知道要加载哪些模块，而不是靠隐式行为碰运气

## 5.3 Phase 3：GLYCAM 扩展模板层补齐

### 目标

建立并逐步缩小 `Xponge` 与 `AmberTools26 GLYCAM_06j` 在扩展糖模板层的差距。

### 文件范围

- [Xponge/forcefield/amber/glycam_06j/terminal.mol2](/mnt/data8t/Software/Xponge/Xponge-origin/Xponge/forcefield/amber/glycam_06j/terminal.mol2)
- `Xponge/forcefield/amber/glycam_06j/*.mol2`
- 可能新增审计脚本与测试文件

### 子阶段划分

#### Phase 3A：覆盖审计脚本

新增一个 `Amber GLYCAM -> Xponge GLYCAM` 的覆盖审计脚本，输入：

- `GLYCAM_06j-1.prep`
- `GLYCAM_amino*.lib`
- `Xponge glycam_06j/*.mol2`

输出四类结果：

1. 已覆盖
2. 名字不同但功能等价
3. 功能团模板缺失
4. 完整修饰糖模板缺失

#### Phase 3B：高优先级功能团与封端块

优先补：

- `SO3`
- `MEX`
- `TBT`
- 必要时 `CA2`

原因：

- 这些是小块模板，结构清晰
- 补齐后对糖修饰/端基支持提升直接
- 实现复杂度明显低于整批特殊糖模板

#### Phase 3C：高价值完整修饰糖模板

第二批补：

- `0TV` 系
- `0AE` 系
- `0AF` 系
- `0GL` / `AGL` 系

这类模板应视为完整 `modified monosaccharide template`，不要误实现成“基础糖 + 名字后缀规则”。

#### Phase 3D：长期补全

当 Phase 3A 到 3C 跑通后，再决定是否追求更高覆盖率。  
这一步建议做成：

- 批量迁移
- 自动审计
- 分族测试

而不是继续手工零散补丁。

### 验收标准

- 有稳定的覆盖报告
- `terminal/function-group` 缺项数量明显下降
- 典型扩展糖模板至少有一批真实可加载例子

## 6. 实现顺序建议

推荐顺序如下：

1. `ff14SB` 的 `HYP -> CHYP` 派生修复
2. `ff19SB` 的 `NHYP` 模板补入
3. `ff19SB` 的 `HYP -> NHYP / CHYP` 派生修复
4. `glycoprotein` 核心模板与说明收口
5. `GLYCAM` 覆盖审计脚本
6. `SO3/MEX/TBT/CA2` 等小块模板
7. `TV/AE/AF/GL` 等完整修饰糖模板

这样可以保证：

- 先解决标准蛋白与糖蛋白核心工作流
- 再推进大规模 GLYCAM 模板扩展

## 7. 测试策略

### 7.1 标准蛋白模板测试

至少新增或补齐：

- `ff14SB` 下 `HYP`
- `ff14SB` 下 `CHYP`
- `ff19SB` 下 `HYP`
- `ff19SB` 下 `NHYP`
- `ff19SB` 下 `CHYP`

### 7.2 糖蛋白核心测试

至少覆盖：

- `NLN`
- `OLS`
- `OLT`
- `OLP`
- `HYP/NHYP/CHYP` 相关组合

### 7.3 扩展模板测试

按族划分：

- 功能团族：
  - `SO3`
  - `MEX`
  - `TBT`
- 特殊糖族：
  - `TV`
  - `AE`
  - `AF`
  - `GL`

### 7.4 审计测试

覆盖审计脚本应进入回归流程，至少保证：

- 不会静默丢失已覆盖模板
- 新增模板会反映在覆盖报告中

## 8. 风险与注意事项

### 8.1 `HYP` 不能机械复用普通残基逻辑

如果直接把 `HYP` 塞进现有 `residues` 列表：

- `ff14SB` 会错误要求 `NHYP`
- 头部几何也会走错分支

### 8.2 `NHYP` 不应手工猜模板

`ff19SB` 的 `NHYP` 应以 Amber 模板为准导入，不建议人工拼接电荷和原子顺序。

### 8.3 扩展糖模板不应过度规则化

像 `0TV`、`0AE`、`0AF`、`0GL` 这类模板更像 Amber 预枚举的完整糖残基，不适合简单抽象成：

- “基础糖 + 取代基后缀”

### 8.4 PDB 自动识别仍是后续问题

即使模板补齐，`load_pdb()` 对普通 PDB 糖命名的自动映射仍然不是本计划的主要目标。  
本计划优先解决“模板是否存在、工作流是否对齐”。

## 9. 工作量评估

- Phase 1：小
- Phase 2：小到中
- Phase 3A：小
- Phase 3B：中
- Phase 3C / 3D：中到大

因此整体判断是：

- 把 `HYP` 和核心糖蛋白修到靠谱，工作量可控
- 把 `Xponge` 推到接近 `AmberTools26 GLYCAM_06j` 完整覆盖，是明确的大工程

## 10. 完成定义

满足以下条件后，可认为该计划的核心目标完成：

1. `ff14SB` 与 `ff19SB` 的 `HYP` 相关行为与 Amber 对齐
2. `ff19SB` 中 `NHYP` 可正常使用
3. `glycoprotein` 核心模板路径稳定
4. `GLYCAM` 扩展模板已有覆盖审计工具
5. 第一批高优先级扩展模板得到补齐

其后，剩余特殊糖模板可转入持续迭代，而不再属于“核心阻塞项”。
