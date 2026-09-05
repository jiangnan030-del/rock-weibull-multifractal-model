# Rock strength Weibull model from multi-fractal damage

## 中文说明

这是针对当前论文《基于损伤多重分形特征的岩石强度 Weibull 参数研究》的一个 Python 复现项目。

目前仓库已经实现：

- Mohr-Coulomb 微元强度公式（文中式 8）
- Weibull 损伤变量（式 7）
- 引入分形参数 `lambda` 的损伤本构关系（式 9）
- 基于峰值条件反求 Weibull 参数 `m` 和 `F0`（式 18-19）
- 根据不同 `lambda` 生成应力-应变曲线
- 整理文中 Table 3 与 Table 4 数据
- 输出 CSV 曲线文件和 matplotlib 图像
- 提供一个可直接运行的 Jupyter Notebook 示例

> 说明：论文页面提供了核心公式和表格，但没有给出所有图的完整实验参数，因此仓库中保留了一个可替换的 demo 参数组，方便后续替换为你自己的实验数据。

### 快速运行

```bash
python3 rock_weibull_model.py
```

这会在仓库目录下生成：

- `data/`：论文表格数据
- `figures/`：示意图
- `curves/`：各组参数对应的应力-应变曲线 CSV

### 命令行用法

#### 1. 只求解给定 `lambda` 对应的 `m` 和 `F0`

```bash
python3 rock_weibull_model.py --mode solve --lambda-values 0.9293 0.9906
```

#### 2. 用自定义峰值参数求解

```bash
python3 rock_weibull_model.py \
  --mode solve \
  --lambda-values 0.9293 \
  --E 15000 \
  --nu 0.25 \
  --sigma3 0 \
  --sigma1c 35 \
  --epsilon1c 0.004 \
  --phi-deg 35
```

#### 3. 导出曲线 CSV

```bash
python3 rock_weibull_model.py --mode curve --lambda-values 0.9293 0.9906
```

### Notebook

- `notebooks/rock_weibull_demo.ipynb`

Notebook 中演示了：

- 如何构造 `PeakState`
- 如何求解 `m` 与 `F0`
- 如何读取论文表格数据
- 如何生成和绘制应力-应变曲线

---

## English summary

This repository is a Python implementation of the paper:

**Liu Shuxin et al. (2011)** — *Weibull distribution parameters of rock strength based on multi-fractal characteristics of rock damage*.

### Implemented features

- micro-element strength based on the Mohr-Coulomb criterion (Eq. 8)
- Weibull damage variable (Eq. 7)
- constitutive relation with fractal damage parameter `lambda` (Eq. 9)
- closed-form solution of Weibull parameters `m` and `F0` from the peak-point conditions (Eqs. 18-19)
- stress-strain curve generation for different crack-spacing / multifractal cases
- digitized Table 3 and Table 4 data
- a Jupyter notebook demo
- curve CSV export and example figures

### Files

- `rock_weibull_model.py`: main implementation
- `data/table3_relative_dimension.csv`: digitized Table 3 data
- `data/table4_singularity.csv`: digitized Table 4 data
- `notebooks/rock_weibull_demo.ipynb`: interactive demo
- `curves/`: generated stress-strain curve CSV files
- `figures/`: generated example plots

### Dependencies

- Python 3.10+
- matplotlib

### Next possible extensions

- fit model parameters against real stress-strain data
- add inverse estimation from experimental curves
- implement multifractal-spectrum calculation directly from crack-network data
- compare relative-dimension and singularity-based `lambda` selection strategies
