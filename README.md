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

## 论文公式汇总

> 记号说明：原文交替使用 $v$ 与 $\nu$，二者均表示泊松比。以下编号按页面中的推导顺序整理；GitHub 支持的 LaTeX 公式可直接渲染。

### 1. 多重分形谱

对盒尺度 $r$，第 $i$ 个盒子的概率为 $p_i(r)$，阶数为 $q$。概率测度族为：

$$
\mu_i(q,r)=\frac{[p_i(r)]^q}{\sum_i[p_i(r)]^q}
$$

支集的 Hausdorff 维数为：

$$
f(q)=\lim_{r\to0}\frac{\sum_i\mu_i(q,r)\lg[\mu_i(q,r)]}{\lg r}
$$

由标度关系：

$$
\alpha_i=\frac{\lg p_i(\delta)}{\lg\delta}
$$

可得奇异性指数：

$$
\alpha(q)=\lim_{r\to0}\frac{\sum_i\mu_i(q,r)\lg[p_i(r)]}{\lg r}
$$

工程计算中，对多个 $r$ 计算上述两个分子，并在无标度区内分别对 $\lg r$ 做线性回归；论文采用相关系数大于 $0.95$ 作为无标度区判据，斜率分别给出 $f(q)$ 与 $\alpha(q)$。

容量维相对分形参数定义为：

$$
\lambda=\frac{D_0}{3},\qquad 0<\lambda<1
$$

其中 $D_0$ 为三维裂纹网络的容量维。论文也尝试直接以多重分形谱奇异性范围（谱宽）作为 $\lambda$。

### 2. 损伤本构关系与 Weibull 分布

Lemaitre 应变等价关系：

$$
[\sigma']=\frac{[\sigma]}{1-D}=\frac{[C][\varepsilon]}{1-D}
$$

面积损伤变量：

$$
D=\frac{A'}{A}
$$

引入分形损伤参数后的关系：

$$
[\sigma']=\frac{[\sigma]}{1-\lambda D}=\frac{[C][\varepsilon]}{1-\lambda D}
$$

岩石微元强度 $F$ 的 Weibull 概率密度：

$$
p(F)=\frac{m}{F_0}\left(\frac{F}{F_0}\right)^{m-1}
\exp\left[-\left(\frac{F}{F_0}\right)^m\right]
$$

由累计破坏概率得到损伤变量：

$$
D=\int_0^F p(x)\,\mathrm{d}x
=1-\exp\left[-\left(\frac{F}{F_0}\right)^m\right]
$$

等围压条件下，Mohr–Coulomb 微元强度为：

$$
F=\frac{E\varepsilon_1\left[(\sigma_1-\sigma_3)-(\sigma_1+\sigma_3)\sin\varphi\right]}
{\sigma_1-2\nu\sigma_3}
$$

对应的轴向应力本构式：

$$
\sigma_1=E\varepsilon_1\left\{\lambda
\exp\left[-\left(\frac{F}{F_0}\right)^m\right]+1-\lambda\right\}
+2\nu\sigma_3
$$

为便于代码生成参数曲线，令

$$
A(F)=\lambda\exp\left[-\left(\frac{F}{F_0}\right)^m\right]+1-\lambda
$$

联立前两式可写成参数形式：

$$
\varepsilon_1(F)=\frac{F-\dfrac{\sigma_3\left[2\nu(1-\sin\varphi)-(1+\sin\varphi)\right]}{A(F)}}
{E(1-\sin\varphi)}
$$

$$
\sigma_1(F)=E\varepsilon_1(F)A(F)+2\nu\sigma_3
$$

### 3. 峰值条件及 Weibull 参数反演

对本构式求导：

$$
\begin{aligned}
\frac{\partial\sigma_1}{\partial\varepsilon_1}
={}&E\left\{\lambda e^{-(F/F_0)^m}+1-\lambda\right\}\\
&+E\varepsilon_1\lambda e^{-(F/F_0)^m}
\left(-\frac{mF^{m-1}}{F_0^m}\right)\\
&\times\left\{
\frac{E[(\sigma_1-\sigma_3)-(\sigma_1+\sigma_3)\sin\varphi]}
{\sigma_1-2\nu\sigma_3}
+\frac{E\varepsilon_1(1-\sin\varphi)}{\sigma_1-2\nu\sigma_3}
\frac{\partial\sigma_1}{\partial\varepsilon_1}
\right.\\
&\qquad\left.
-\frac{E\varepsilon_1[(\sigma_1-\sigma_3)-(\sigma_1+\sigma_3)\sin\varphi]}
{(\sigma_1-2\nu\sigma_3)^2}
\frac{\partial\sigma_1}{\partial\varepsilon_1}
\right\}.
\end{aligned}
$$

峰值点满足：

$$
\sigma_1=\sigma_{1c},\qquad
\varepsilon_1=\varepsilon_{1c},\qquad
F=F_c,\qquad
\left.\frac{\partial\sigma_1}{\partial\varepsilon_1}\right|_c=0
$$

因此：

$$
\begin{aligned}
0={}&E\left\{\lambda e^{-(F_c/F_0)^m}+1-\lambda\right\}\\
&+E\varepsilon_{1c}\lambda e^{-(F_c/F_0)^m}
\left(-\frac{mF_c^{m-1}}{F_0^m}\right)
\frac{E[(\sigma_{1c}-\sigma_3)-(\sigma_{1c}+\sigma_3)\sin\varphi]}
{\sigma_{1c}-2\nu\sigma_3}.
\end{aligned}
$$

峰值微元强度：

$$
F_c=\frac{E\varepsilon_{1c}
[(\sigma_{1c}-\sigma_3)-(\sigma_{1c}+\sigma_3)\sin\varphi]}
{\sigma_{1c}-2\nu\sigma_3}
$$

代入后得到：

$$
\lambda e^{-(F_c/F_0)^m}+1-\lambda
-\lambda e^{-(F_c/F_0)^m}\frac{mF_c^m}{F_0^m}=0
$$

峰值本构关系：

$$
\sigma_{1c}=E\varepsilon_{1c}
\left\{\lambda e^{-(F_c/F_0)^m}+1-\lambda\right\}
+2\nu\sigma_3
$$

从而有：

$$
\lambda e^{-(F_c/F_0)^m}+1-\lambda
=\frac{\sigma_{1c}-2\nu\sigma_3}{E\varepsilon_{1c}}
$$

$$
\lambda e^{-(F_c/F_0)^m}
=\frac{\sigma_{1c}-2\nu\sigma_3}{E\varepsilon_{1c}}+\lambda-1
$$

$$
\frac{mF_c^m}{F_0^m}
=\frac{\sigma_{1c}-2\nu\sigma_3}
{\sigma_{1c}-2\nu\sigma_3+(\lambda-1)E\varepsilon_{1c}}
$$

令

$$
R=\frac{\sigma_{1c}-2\nu\sigma_3+(\lambda-1)E\varepsilon_{1c}}
{E\varepsilon_{1c}\lambda}
$$

则：

$$
\left(\frac{F_c}{F_0}\right)^m=-\ln R
$$

$$
m=-\frac{\sigma_{1c}-2\nu\sigma_3}
{\sigma_{1c}-2\nu\sigma_3+(\lambda-1)E\varepsilon_{1c}}
\frac{1}{\ln R}
$$

页面未单独写出但代码由上一式直接使用的 $F_0$ 显式解为：

$$
F_0=\frac{F_c}{[-\ln R]^{1/m}}
$$

物理解要求：

$$
0<\lambda<1,\qquad 0<R<1,\qquad m>0,\qquad F_0>0
$$

## 项目审查结果与复现边界

### 已覆盖

- Weibull 密度、损伤累计概率、Mohr–Coulomb 微元强度和分形损伤本构式。
- 峰值点 $F_c$、$m$、$F_0$ 的反演，以及参数化应力–应变曲线。
- 论文表 3、表 4 的数据、曲线 CSV、命令行入口和 Notebook 示例。

### 本次补正

- README 已补齐本页的多重分形、本构、峰值条件和参数反演公式。
- 增加表 2 的容量维与奇异性范围数据文件。
- 参数求解增加 $0<\lambda<1$、$0<R<1$ 及正参数校验，避免产生非物理解或复数结果。

### 尚未实现或无法严格复现

- 表 1 的原始裂纹产状数据在页面 OCR 中错位，不能可靠还原为机器可读表。
- 尚未实现论文第 1 节的三维裂纹 Monte Carlo 网络生成。
- 尚未实现从盒计数数据自动选择无标度区并计算 $f(q)$、$\alpha(q)$ 的完整流水线；README 现已保留公式和判据。
- 图 6、图 7 所述的单独 $m$ / $F_0$ 敏感性扫描没有独立 CLI 模式。
- 论文页面未提供完整原始应力–应变实验数据和全部试验参数，因此仓库默认参数生成的是说明性曲线，不应视为论文原图的逐点复现。
- 表 3/表 4 的已发表 $m,F_0$ 与默认 demo 峰值参数来自不同信息层；混合使用仅用于展示模型行为。严格验证应输入论文对应试样的 $E,\nu,\sigma_3,\sigma_{1c},\varepsilon_{1c},\varphi$。
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

### Example figures

SVG exports are included for easier preview on GitHub:

- `figures/table3_m_vs_spacing.svg`
- `figures/table3_F0_vs_spacing.svg`
- `figures/table4_m_vs_spacing.svg`
- `figures/table4_F0_vs_spacing.svg`
- `figures/stress_strain_table3.svg`
- `figures/stress_strain_table4.svg`

### Dependencies

- Python 3.10+
- matplotlib
- pandas

### Next possible extensions

- fit model parameters against real stress-strain data
- add inverse estimation from experimental curves
- implement multifractal-spectrum calculation directly from crack-network data
- compare relative-dimension and singularity-based `lambda` selection strategies
