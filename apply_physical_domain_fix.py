from pathlib import Path

readme = Path('README.md')
s = readme.read_text(encoding='utf-8')
s = s.replace('\\lambda=\\frac{D_0}{3},\\qquad 0<\\lambda<1', '\\lambda=\\frac{D_0}{3},\\qquad 0\\lt\\lambda\\le 1')
old = '''物理解要求：

$$
0<\\lambda<1,\\qquad 0<R<1,\\qquad m>0,\\qquad F_0>0
$$

> 原文交替使用 $v$ 和 $\\nu$，本仓库统一使用 `nu` 表示泊松比。'''
new = '''为明确物理解的输入范围，定义：

$$
A=\\sigma_{1c}-2\\nu\\sigma_3,\\qquad
B=A+(\\lambda-1)E\\varepsilon_{1c}
$$

$$
Q=(\\sigma_{1c}-\\sigma_3)-(\\sigma_{1c}+\\sigma_3)\\sin\\varphi
$$

物理解要求：

$$
E\\gt0,\\quad \\varepsilon_{1c}\\gt0,\\quad 0\\lt\\lambda\\le1
$$

$$
0\\le\\nu\\lt0.5,\\quad 0\\le\\varphi\\lt90^\\circ
$$

$$
A\\gt0,\\quad B\\gt0,\\quad Q\\gt0,\\quad 0\\lt R\\lt1
$$

在这些输入条件下，$F_c\\gt0$，求解结果应满足：

$$
m\\gt0,\\qquad F_0\\gt0
$$

由于 $E\\varepsilon_{1c}\\lambda\\gt0$，$0\\lt R\\lt1$ 等价于：

$$
(1-\\lambda)E\\varepsilon_{1c}
\\lt \\sigma_{1c}-2\\nu\\sigma_3
\\lt E\\varepsilon_{1c}
$$

> 论文正文采用 $0\\lt\\lambda\\lt1$。从 $\\lambda=D_0/3$ 和反演公式的定义域看，$\\lambda=1$ 可作为 $D_0=3$ 的极限情形；本仓库允许该边界值，但仍要求 $0\\lt R\\lt1$。原文交替使用 $v$ 和 $\\nu$，本仓库统一使用 `nu`。'''
if old not in s:
    raise SystemExit('README physical-domain block not found')
s = s.replace(old, new, 1)
readme.write_text(s, encoding='utf-8')

code_path = Path('rock_weibull_model.py')
code = code_path.read_text(encoding='utf-8')
code = code.replace('A physical Weibull solution requires 0 < lambda < 1 and 0 < R < 1,', 'A physical Weibull solution requires 0 < lambda <= 1 and 0 < R < 1,')
old_code = '''    if not 0.0 < lambda_value < 1.0:
        raise ValueError("lambda must satisfy 0 < lambda < 1")

    numerator = peak.sigma1c - 2.0 * peak.nu * peak.sigma3
'''
new_code = '''    if peak.E <= 0.0 or peak.epsilon1c <= 0.0:
        raise ValueError("E and epsilon1c must be positive")
    if not 0.0 <= peak.nu < 0.5:
        raise ValueError("nu must satisfy 0 <= nu < 0.5")
    if not 0.0 <= peak.phi_deg < 90.0:
        raise ValueError("phi_deg must satisfy 0 <= phi_deg < 90")
    if not 0.0 < lambda_value <= 1.0:
        raise ValueError("lambda must satisfy 0 < lambda <= 1")

    numerator = peak.sigma1c - 2.0 * peak.nu * peak.sigma3
    if numerator <= 0.0:
        raise ValueError("sigma1c - 2*nu*sigma3 must be positive")
    phi = math.radians(peak.phi_deg)
    strength_factor = (peak.sigma1c - peak.sigma3) - (peak.sigma1c + peak.sigma3) * math.sin(phi)
    if strength_factor <= 0.0:
        raise ValueError("peak Mohr-Coulomb strength factor must be positive")
'''
if old_code not in code:
    raise SystemExit('solver validation block not found')
code_path.write_text(code.replace(old_code, new_code, 1), encoding='utf-8')
