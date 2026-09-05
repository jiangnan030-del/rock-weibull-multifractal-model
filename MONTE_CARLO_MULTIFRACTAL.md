# 三维微裂纹 Monte Carlo 与多重分形计算

`microcrack_multifractal.py` 完成论文第 1～2 节的可执行流程：

1. 按校正后的表 1 生成两组裂纹。
2. 裂纹中心在立方体中均匀分布。
3. 倾向和倾角按正态分布抽样，标准差取表中方差的平方根。
4. 迹长按负指数分布抽样；每条裂纹表示为厚度忽略的圆盘。
5. 在圆盘面积上均匀采样点并裁剪到模拟立方体。
6. 对多个三维盒尺度统计占据盒及概率 $p_i(r)$。
7. 计算 Chhabra–Jensen 测度 $\mu_i(q,r)$。
8. 自动搜索连续无标度区；要求 $\alpha$ 与 $f$ 两个回归的相关系数绝对值均不低于阈值，默认 0.95。
9. 输出 $\alpha(q)$、$f(q)$、容量维 $D_0=f(0)$、相对维数 $\lambda=D_0/3$ 和谱宽。

## 快速运行

```bash
python3 microcrack_multifractal.py
```

```bash
python3 microcrack_multifractal.py --spacing 0.5 --cube-size 10 \
  --points-per-crack 300 --divisions 4,5,6,8,10,12,16,20 \
  --q-min -5 --q-max 5 --q-step 0.5 --correlation-threshold 0.95
```

批量计算论文五个间距：

```bash
for spacing in 0.1 0.3 0.5 0.7 0.9; do
  python3 microcrack_multifractal.py --spacing "$spacing"
done
```

## 输出

默认写入 `outputs/multifractal_spacing_<间距>/`：

- `microcrack_points.csv`：裂纹面采样点、编号和裂纹组。
- `box_count_scales.csv`：盒尺度与占据盒数。
- `multifractal_spectrum.csv`：每个 $q$ 的 $\alpha(q)$、$f(q)$、自动选择的无标度区、相关系数和 $R^2$。
- `summary.json`：$D_0$、$\lambda$、谱宽、随机种子和参数。
- `microcrack_network_3d.png`：三维微裂纹网络。
- `multifractal_spectrum.png`：$f(\alpha)$ 谱。

## 自动无标度区

程序枚举不少于 `--min-window-points` 的连续尺度窗口，同时拟合
$\sum_i\mu_i\ln p_i$ 与 $\sum_i\mu_i\ln\mu_i$ 关于 $\ln r$ 的直线。优先选择两个相关系数绝对值都不低于 0.95、斜率满足三维物理范围的最长窗口；若没有窗口通过，返回最佳候选并将 `threshold_passed=false`，不静默伪装为有效结果。

## 建模边界

论文未完整说明裂纹数量、圆盘离散密度、边界截断方式和负指数分布参数化。本实现把这些选择显式参数化并记录随机种子。表 1 的迹长“均值、方差、负指数分布”并不满足标准指数分布方差必须等于均值平方的约束，因此实现使用表中均值作为指数尺度，同时保留原方差，不虚构额外参数。
