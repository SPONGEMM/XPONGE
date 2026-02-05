# Xponge 非正交盒角度支持改动说明

本文档总结了为支持非正交盒角度读写与构建所做的代码修改与使用方式。

## 目标
- 读取文件时保留盒子角度信息（PDB/GRO/CIF）。
- 保存 PDB 时输出 CRYST1 角度信息。
- 构建晶体时允许非正交晶胞（角度不再限制 <=90）。
- 盒子角度是否读入为可选项。

## 主要改动

### 1) 盒角度解析工具
**文件**：`Xponge/helper/math.py`
- 新增 `get_length_angle_from_basis_vectors(v1, v2, v3)`：
  - 由三基矢计算 `[a,b,c]` 与 `[alpha,beta,gamma]`。

并在 `Xponge/helper/__init__.py` 中导出该函数。

---

### 2) GRO 读取支持 9 值 triclinic box
**文件**：`Xponge/load.py`
- `load_gro` 新增参数 `read_box_angle=True`。
- 支持 9 值 box（triclinic），计算长度/角度并写入 `mol.box_length / mol.box_angle`。
- 若 `read_box_angle=False`，角度强制为 `[90,90,90]`。

**新签名**：
```python
load_gro(filename, mol=None, read_box_angle=True)
```

---

### 3) PDB 读取 CRYST1
**文件**：`Xponge/load.py`
- `load_pdb` 新增参数 `read_cryst1=True`。
- 若读到 `CRYST1`，设置 `molecule.box_length` 与 `molecule.box_angle`。

**新签名**：
```python
load_pdb(file, ..., ignore_conect=True, read_cryst1=True)
```

---

### 4) CIF 读取角度可选
**文件**：`Xponge/assign/__init__.py`
- `get_assignment_from_cif` 新增参数 `keep_cell_angle=True`。
- `orthogonal_threshold` 默认改为 `None`，表示不再自动“贴近 90°”。

**新签名**：
```python
get_assignment_from_cif(file, total_charge=0, orthogonal_threshold=None, keep_cell_angle=True)
```

---

### 5) PDB 保存 CRYST1
**文件**：`Xponge/build.py`
- `save_pdb` 新增参数 `write_cryst1=True`。
- 若 `mol.box_length` 存在，则输出 `CRYST1` 并写入角度。

**新签名**：
```python
save_pdb(cls, filename=None, write_cryst1=True)
```

---

### 6) SPONGE 输入保留角度
**文件**：`Xponge/__init__.py`
- `_do_initial` 将 `pbc_box` 从 3 长度扩展为 6 值 `[a,b,c,alpha,beta,gamma]`。
- 如果未设置 `box_angle`，默认 `[90,90,90]`。

---

### 7) 非正交晶胞构建
**文件**：`Xponge/process.py`
- `Lattice.__init__` 取消 `angle <= 90` 限制，改为 `(0,180)`。
- `Lattice.create` 允许 `PrismRegion` 作为 box，支持非正交盒构建。

---

## 使用示例

### 读取 PDB 并保留角度
```python
pro = Xponge.load_pdb("your.pdb", read_cryst1=True)
```

### 读取 GRO triclinic box
```python
crd, box = Xponge.load_gro("your.gro", read_box_angle=True)
```

### CIF 保留角度
```python
assign, lattice_info = Xponge.get_assignment_from_cif("your.cif", keep_cell_angle=True)
```

### 保存 PDB 输出 CRYST1
```python
Xponge.save_pdb(pro, "out.pdb", write_cryst1=True)
```

### 非正交晶格构建
```python
from Xponge.process import Lattice, PrismRegion

box = PrismRegion(0,0,0, vx,vy,vz, wx,wy,wz, ux,uy,uz, boundary=True)
lat = Lattice(basis_molecule=restype, cell_length=[a,b,c], cell_angle=[alpha,beta,gamma])
mol = lat.create(box, box)
```

---

## 兼容性说明
- 正交盒流程保持兼容。
- 若外部引擎仅接受 3 值盒长，需自行确认其支持 6 值 `[a,b,c,alpha,beta,gamma]`。

## 可选后续增强
- `save_gro` 支持输出 9 值 triclinic 盒。
- `save_sponge_input` 明确写入角度字段的规范化格式说明。

