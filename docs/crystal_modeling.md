# Xponge 晶体建模指南（含苯晶体示例）

本文介绍如何用 Xponge 构建正交与非正交晶体，并给出一个**苯晶体**的完整示例：
- SMILES 构建小分子
- GAFF 力场赋型/电荷
- 非正交晶胞晶体构建

> 依赖：SMILES 解析需要 RDKit（Xponge 会调用 RDKit）。

---

## 一、通用流程

1. **准备基元分子（basis_molecule）**
   - 通过 SMILES、mol2 或 CIF 得到 `Assign` 或 `ResidueType`。
   - 用力场规则（如 GAFF）分配原子类型和电荷。
   - 转为 `ResidueType` 或 `Molecule` 作为晶胞基元。

2. **设置晶胞参数**
   - 正交：`cell_angle=[90,90,90]`
   - 非正交：指定实际角度 `[alpha, beta, gamma]`

3. **定义晶体区域**
   - 正交盒：`BlockRegion`
   - 非正交盒：`PrismRegion`（定义三条盒矢）

4. **构建晶体**
   - 用 `Lattice.create(box, region)` 进行复制填充。

---

## 二、Lattice 相关功能说明

`Lattice` 用来定义晶格复制规则（晶胞参数 + 基元 + 晶胞内位置）。

**核心参数：**
- `style`：晶格模板。`"custom"` 表示完全自定义；`"template:NAME"` 可注册模板样式；也可用内置 `"sc" "fcc" "bcc" "hcp" "diamond"`。
- `basis_molecule`：晶胞基元，可为 `ResidueType / Residue / Molecule`，坐标即基元内部坐标。
- `cell_length`：晶胞长度 `[a,b,c]`（Å）。
- `cell_angle`：晶胞角 `[alpha,beta,gamma]`（度）。
- `basis_position`：晶胞内基元位置列表，每个位置是 **晶胞基矢坐标**（不是笛卡尔坐标）。
- `scale`：**非模板 lattice 必填**，用于整体缩放晶胞基矢与基元（常用 `1.0`）。
- `spacing`：在基矢方向额外留白，用于稀疏填充。

**常见坑：**
- `style="custom"` 时必须提供 `scale`，否则会报：
  `ValueError: scale should not be None for a non-template lattice`

---

## 二、苯晶体（非正交）示例

### 1) 从 SMILES 生成苯
```python
import Xponge
import Xponge.forcefield.amber.gaff as gaff
from Xponge.process import Lattice, PrismRegion
from Xponge.helper.math import get_basis_vectors_from_length_and_angle

# 1. SMILES 构建 Assign
assign = Xponge.get_assignment_from_smiles("c1ccccc1")

# 2. GAFF 赋原子类型
assign.Determine_Atom_Type("gaff")

# 3. 计算电荷（可选：若已有电荷可跳过）
assign.Calculate_Charge("RESP")

# 4. 转 ResidueType 作为晶胞基元
ben = assign.toResidueType("BEN")
```

### 2) 设置非正交晶胞
假设晶胞参数（示例）为：
- `a=7.5, b=7.5, c=5.0`
- `alpha=90, beta=90, gamma=60`（非正交）

```python
# 通过长度/角度生成盒矢
basis = get_basis_vectors_from_length_and_angle(7.5, 7.5, 5.0, 90, 90, 60)
# basis 是 3x3，每一行是一条盒矢
v1 = basis[0]
v2 = basis[1]
v3 = basis[2]
```

### 3) 构建非正交晶体
```python
# 使用 PrismRegion 定义非正交盒
box = PrismRegion(
    0, 0, 0,
    v1[0], v1[1], v1[2],
    v2[0], v2[1], v2[2],
    v3[0], v3[1], v3[2],
    boundary=False  # 避免边界重复引入一排
)

# Lattice 构建
lattice = Lattice(
    style="custom",
    basis_molecule=ben,
    scale=1.0,  # 非模板 lattice 必填
    cell_length=[7.5, 7.5, 5.0],
    cell_angle=[90, 90, 60],
    basis_position=[[0, 0, 0]]
)

crystal = lattice.create(box, box)

# 保存结果
Xponge.save_pdb(crystal, "benzene_crystal.pdb")
```

---

## 三、说明与注意事项

- **基元坐标来源**：`basis_molecule` 本身的原子坐标。
- **basis_position**：晶胞内多个基元的相对位置（以晶胞基矢坐标系表示）。
- **非正交盒角度**：通过 `PrismRegion` + `cell_angle` 保留。
- **电荷**：小分子建议用 RESP；若已有电荷可跳过。

---

## 四、正交晶体对照
若晶胞为正交（90°）：
```python
from Xponge.process import BlockRegion

box = BlockRegion(0, 0, 0, 30, 30, 30, boundary=True)

lattice = Lattice(
    style="custom",
    basis_molecule=ben,
    scale=1.0,  # 非模板 lattice 必填
    cell_length=[7.5, 7.5, 5.0],
    cell_angle=[90, 90, 90],
    basis_position=[[0, 0, 0]]
)

crystal = lattice.create(box, box)
```

---

## 五、常见问题

**Q1: 为什么要用 PrismRegion？**
> BlockRegion 仅支持轴对齐正交盒；非正交必须用 PrismRegion。

**Q2: 基元要用 Molecule 还是 ResidueType？**
> 两者都可以。ResidueType 更轻量，Molecule 可保留更完整结构。

**Q3: RDKit 未安装怎么办？**
> SMILES 路径依赖 RDKit。可以改用 mol2/CIF 作为输入。
