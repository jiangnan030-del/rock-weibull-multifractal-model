from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Iterable
import csv
import math

import matplotlib.pyplot as plt


@dataclass
class PeakState:
    """Peak-point parameters from the constitutive derivation."""

    E: float  # Elastic modulus, MPa
    nu: float  # Poisson ratio
    sigma3: float  # Confining pressure, MPa
    sigma1c: float  # Peak axial stress, MPa
    epsilon1c: float  # Peak axial strain
    phi_deg: float  # Friction angle, degrees


@dataclass
class WeibullParameters:
    m: float
    F0: float


@dataclass
class RelativeDimensionRow:
    spacing_cm: float
    lambda_value: float
    m: float
    F0_MPa: float


TABLE3_RELATIVE_DIMENSION: list[RelativeDimensionRow] = [
    RelativeDimensionRow(0.9, 0.9293, 2.2366, 32.0727),
    RelativeDimensionRow(0.7, 0.9364, 2.2307, 32.2401),
    RelativeDimensionRow(0.5, 0.9533, 2.2170, 32.6409),
    RelativeDimensionRow(0.3, 0.9906, 2.1893, 33.5194),
    RelativeDimensionRow(0.1, 0.9991, 2.1833, 33.7199),
]

TABLE4_SINGULARITY: list[RelativeDimensionRow] = [
    RelativeDimensionRow(0.9, 0.6578, 2.6634, 25.4432),
    RelativeDimensionRow(0.7, 0.5130, 3.4497, 22.0651),
    RelativeDimensionRow(0.5, 0.5662, 3.0397, 23.2395),
    RelativeDimensionRow(0.3, 0.8165, 2.3544, 29.3504),
    RelativeDimensionRow(0.1, 0.9113, 2.2525, 31.6424),
]


def mohr_coulomb_microelement_strength(
    epsilon1: float,
    sigma1: float,
    sigma3: float,
    E: float,
    nu: float,
    phi_deg: float,
) -> float:
    """Equation (8) in the paper."""
    phi = math.radians(phi_deg)
    numerator = E * epsilon1 * ((sigma1 - sigma3) - (sigma1 + sigma3) * math.sin(phi))
    denominator = sigma1 - 2.0 * nu * sigma3
    if abs(denominator) < 1e-12:
        raise ZeroDivisionError("sigma1 - 2*nu*sigma3 is too close to zero")
    return numerator / denominator



def solve_weibull_parameters(peak: PeakState, lambda_value: float) -> WeibullParameters:
    """Solve m and F0 from equations (19) and (18)."""
    numerator = peak.sigma1c - 2.0 * peak.nu * peak.sigma3
    denom = numerator + (lambda_value - 1.0) * peak.E * peak.epsilon1c
    if abs(denom) < 1e-12:
        raise ZeroDivisionError("Equation denominator is too close to zero")

    log_argument = denom / (peak.E * peak.epsilon1c * lambda_value)
    if log_argument <= 0.0 or math.isclose(log_argument, 1.0, rel_tol=0.0, abs_tol=1e-12):
        raise ValueError(
            "Invalid logarithm argument from the peak-state inputs; choose physically consistent values"
        )

    m = -(numerator / denom) / math.log(log_argument)

    Fc = mohr_coulomb_microelement_strength(
        epsilon1=peak.epsilon1c,
        sigma1=peak.sigma1c,
        sigma3=peak.sigma3,
        E=peak.E,
        nu=peak.nu,
        phi_deg=peak.phi_deg,
    )
    ratio_power_m = -math.log(log_argument)
    F0 = Fc / (ratio_power_m ** (1.0 / m))
    return WeibullParameters(m=m, F0=F0)



def damage_variable(F: float, m: float, F0: float) -> float:
    return 1.0 - math.exp(-((F / F0) ** m))



def axial_stress_from_F(F: float, peak: PeakState, lambda_value: float, params: WeibullParameters) -> tuple[float, float]:
    """Return (epsilon1, sigma1) from equations (8) and (9).

    Rearranging Eq. (8) after substituting Eq. (9) gives:
        F = E*epsilon1*(1-sin(phi)) + sigma3*(2*nu*(1-sin(phi))-(1+sin(phi))) / A
    where A = lambda*exp(-(F/F0)^m) + 1 - lambda.
    """
    phi = math.radians(peak.phi_deg)
    sin_phi = math.sin(phi)
    factor = lambda_value * math.exp(-((F / params.F0) ** params.m)) + 1.0 - lambda_value
    if abs(factor) < 1e-12:
        raise ZeroDivisionError("Damage factor is too close to zero")

    correction = peak.sigma3 * (2.0 * peak.nu * (1.0 - sin_phi) - (1.0 + sin_phi)) / factor
    epsilon1 = (F - correction) / (peak.E * (1.0 - sin_phi))
    sigma1 = peak.E * epsilon1 * factor + 2.0 * peak.nu * peak.sigma3
    return epsilon1, sigma1



def stress_strain_curve(
    peak: PeakState,
    lambda_value: float,
    params: WeibullParameters,
    n: int = 300,
    F_max_factor: float = 2.0,
) -> list[tuple[float, float]]:
    Fc = mohr_coulomb_microelement_strength(
        epsilon1=peak.epsilon1c,
        sigma1=peak.sigma1c,
        sigma3=peak.sigma3,
        E=peak.E,
        nu=peak.nu,
        phi_deg=peak.phi_deg,
    )
    values: list[tuple[float, float]] = []
    last_epsilon = -float("inf")
    for i in range(1, n + 1):
        F = Fc * F_max_factor * i / n
        epsilon1, sigma1 = axial_stress_from_F(F, peak, lambda_value, params)
        if math.isfinite(epsilon1) and math.isfinite(sigma1) and epsilon1 > 0.0:
            if epsilon1 >= last_epsilon:
                values.append((epsilon1, sigma1))
                last_epsilon = epsilon1
    return values



def save_rows_csv(rows: Iterable[RelativeDimensionRow], path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f)
        writer.writerow(["spacing_cm", "lambda", "m", "F0_MPa"])
        for row in rows:
            writer.writerow([row.spacing_cm, row.lambda_value, row.m, row.F0_MPa])



def plot_parameter_vs_spacing(rows: list[RelativeDimensionRow], attribute: str, ylabel: str, out_path: Path) -> None:
    out_path.parent.mkdir(parents=True, exist_ok=True)
    xs = [r.spacing_cm for r in rows]
    ys = [getattr(r, attribute) for r in rows]
    plt.figure(figsize=(6, 4))
    plt.plot(xs, ys, marker="o")
    plt.gca().invert_xaxis()
    plt.xlabel("Micro-crack spacing (cm)")
    plt.ylabel(ylabel)
    plt.grid(True, linestyle="--", alpha=0.4)
    plt.tight_layout()
    plt.savefig(out_path, dpi=200)
    plt.close()



def plot_stress_strain_family(
    peak: PeakState,
    rows: list[RelativeDimensionRow],
    title: str,
    out_path: Path,
) -> None:
    out_path.parent.mkdir(parents=True, exist_ok=True)
    plt.figure(figsize=(6, 4))
    for row in rows:
        params = WeibullParameters(m=row.m, F0=row.F0_MPa)
        curve = stress_strain_curve(peak, row.lambda_value, params)
        if curve:
            xs = [eps for eps, _ in curve]
            ys = [sig for _, sig in curve]
            plt.plot(xs, ys, label=f"spacing={row.spacing_cm} cm")
    plt.xlabel("Axial strain")
    plt.ylabel("Axial stress (MPa)")
    plt.title(title)
    plt.grid(True, linestyle="--", alpha=0.4)
    handles, labels = plt.gca().get_legend_handles_labels()
    if labels:
        plt.legend(fontsize=8)
    plt.tight_layout()
    plt.savefig(out_path, dpi=200)
    plt.close()



def build_demo_outputs(base_dir: Path) -> None:
    data_dir = base_dir / "data"
    fig_dir = base_dir / "figures"
    save_rows_csv(TABLE3_RELATIVE_DIMENSION, data_dir / "table3_relative_dimension.csv")
    save_rows_csv(TABLE4_SINGULARITY, data_dir / "table4_singularity.csv")

    plot_parameter_vs_spacing(TABLE3_RELATIVE_DIMENSION, "m", "m", fig_dir / "table3_m_vs_spacing.png")
    plot_parameter_vs_spacing(TABLE3_RELATIVE_DIMENSION, "F0_MPa", "F0 (MPa)", fig_dir / "table3_F0_vs_spacing.png")
    plot_parameter_vs_spacing(TABLE4_SINGULARITY, "m", "m", fig_dir / "table4_m_vs_spacing.png")
    plot_parameter_vs_spacing(TABLE4_SINGULARITY, "F0_MPa", "F0 (MPa)", fig_dir / "table4_F0_vs_spacing.png")

    # Demo peak-state parameters; replace with experiment-specific values when available.
    demo_peak = PeakState(E=15000.0, nu=0.25, sigma3=0.0, sigma1c=35.0, epsilon1c=0.004, phi_deg=35.0)
    plot_stress_strain_family(
        demo_peak,
        TABLE3_RELATIVE_DIMENSION,
        "Stress-strain family from Table 3",
        fig_dir / "stress_strain_table3.png",
    )
    plot_stress_strain_family(
        demo_peak,
        TABLE4_SINGULARITY,
        "Stress-strain family from Table 4",
        fig_dir / "stress_strain_table4.png",
    )


if __name__ == "__main__":
    project_dir = Path(__file__).resolve().parent
    build_demo_outputs(project_dir)
    demo_peak = PeakState(E=15000.0, nu=0.25, sigma3=0.0, sigma1c=35.0, epsilon1c=0.004, phi_deg=35.0)
    for lambda_value in [0.9293, 0.9906, 0.6578]:
        params = solve_weibull_parameters(demo_peak, lambda_value)
        print(f"lambda={lambda_value:.4f} -> m={params.m:.4f}, F0={params.F0:.4f} MPa")
    print("Demo outputs written to:", project_dir)
