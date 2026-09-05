from __future__ import annotations

import argparse
import csv
import json
import math
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

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


def default_demo_peak() -> PeakState:
    """A configurable demo peak state for reproducing example outputs.

    Replace these values with experiment-specific parameters when available.
    """
    return PeakState(E=15000.0, nu=0.25, sigma3=0.0, sigma1c=35.0, epsilon1c=0.004, phi_deg=35.0)


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
    """Solve m and F0 from equations (19) and (18).

    A physical Weibull solution requires 0 < lambda <= 1 and 0 < R < 1,
    where R is the logarithm argument in equation (18).
    """
    if peak.E <= 0.0 or peak.epsilon1c <= 0.0:
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
    denom = numerator + (lambda_value - 1.0) * peak.E * peak.epsilon1c
    if abs(denom) < 1e-12:
        raise ZeroDivisionError("Equation denominator is too close to zero")

    log_argument = denom / (peak.E * peak.epsilon1c * lambda_value)
    if not 0.0 < log_argument < 1.0:
        raise ValueError(
            "Physical solution requires the equation (18) logarithm argument R to satisfy 0 < R < 1"
        )

    m = -(numerator / denom) / math.log(log_argument)
    if not math.isfinite(m) or m <= 0.0:
        raise ValueError("Solved Weibull shape parameter m must be positive and finite")

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
    if not math.isfinite(F0) or F0 <= 0.0:
        raise ValueError("Solved Weibull scale parameter F0 must be positive and finite")
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


def save_curve_csv(curve: list[tuple[float, float]], path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f)
        writer.writerow(["epsilon1", "sigma1_MPa"])
        for epsilon1, sigma1 in curve:
            writer.writerow([epsilon1, sigma1])


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


def build_demo_outputs(base_dir: Path, peak: PeakState | None = None) -> None:
    peak = peak or default_demo_peak()
    data_dir = base_dir / "data"
    fig_dir = base_dir / "figures"
    curve_dir = base_dir / "curves"

    save_rows_csv(TABLE3_RELATIVE_DIMENSION, data_dir / "table3_relative_dimension.csv")
    save_rows_csv(TABLE4_SINGULARITY, data_dir / "table4_singularity.csv")

    plot_parameter_vs_spacing(TABLE3_RELATIVE_DIMENSION, "m", "m", fig_dir / "table3_m_vs_spacing.png")
    plot_parameter_vs_spacing(TABLE3_RELATIVE_DIMENSION, "F0_MPa", "F0 (MPa)", fig_dir / "table3_F0_vs_spacing.png")
    plot_parameter_vs_spacing(TABLE4_SINGULARITY, "m", "m", fig_dir / "table4_m_vs_spacing.png")
    plot_parameter_vs_spacing(TABLE4_SINGULARITY, "F0_MPa", "F0 (MPa)", fig_dir / "table4_F0_vs_spacing.png")

    plot_stress_strain_family(
        peak,
        TABLE3_RELATIVE_DIMENSION,
        "Stress-strain family from Table 3",
        fig_dir / "stress_strain_table3.png",
    )
    plot_stress_strain_family(
        peak,
        TABLE4_SINGULARITY,
        "Stress-strain family from Table 4",
        fig_dir / "stress_strain_table4.png",
    )

    for dataset_name, rows in {
        "table3": TABLE3_RELATIVE_DIMENSION,
        "table4": TABLE4_SINGULARITY,
    }.items():
        for row in rows:
            params = WeibullParameters(m=row.m, F0=row.F0_MPa)
            curve = stress_strain_curve(peak, row.lambda_value, params)
            file_name = f"{dataset_name}_spacing_{str(row.spacing_cm).replace('.', '_')}.csv"
            save_curve_csv(curve, curve_dir / file_name)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Python implementation of the rock-strength Weibull model based on multifractal damage."
    )
    parser.add_argument("--mode", choices=["demo", "solve", "curve"], default="demo")
    parser.add_argument("--outdir", default=None, help="Output directory. Defaults to the script directory.")
    parser.add_argument("--lambda-values", nargs="+", type=float, default=[0.9293, 0.9906, 0.6578])
    parser.add_argument("--curve-points", type=int, default=300)
    parser.add_argument("--curve-max-factor", type=float, default=2.0)

    parser.add_argument("--E", type=float, default=15000.0)
    parser.add_argument("--nu", type=float, default=0.25)
    parser.add_argument("--sigma3", type=float, default=0.0)
    parser.add_argument("--sigma1c", type=float, default=35.0)
    parser.add_argument("--epsilon1c", type=float, default=0.004)
    parser.add_argument("--phi-deg", type=float, default=35.0)
    return parser.parse_args()


def peak_from_args(args: argparse.Namespace) -> PeakState:
    return PeakState(
        E=args.E,
        nu=args.nu,
        sigma3=args.sigma3,
        sigma1c=args.sigma1c,
        epsilon1c=args.epsilon1c,
        phi_deg=args.phi_deg,
    )


def run_solve_mode(peak: PeakState, lambda_values: list[float]) -> None:
    results = []
    for lambda_value in lambda_values:
        params = solve_weibull_parameters(peak, lambda_value)
        results.append(
            {
                "lambda": lambda_value,
                "m": params.m,
                "F0_MPa": params.F0,
            }
        )
    print(json.dumps(results, ensure_ascii=False, indent=2))


def run_curve_mode(
    peak: PeakState,
    lambda_values: list[float],
    n: int,
    F_max_factor: float,
    outdir: Path,
) -> None:
    outdir.mkdir(parents=True, exist_ok=True)
    for lambda_value in lambda_values:
        params = solve_weibull_parameters(peak, lambda_value)
        curve = stress_strain_curve(peak, lambda_value, params, n=n, F_max_factor=F_max_factor)
        file_name = f"curve_lambda_{str(lambda_value).replace('.', '_')}.csv"
        save_curve_csv(curve, outdir / file_name)
        print(f"Saved {file_name}")


def main() -> None:
    args = parse_args()
    peak = peak_from_args(args)
    project_dir = Path(args.outdir).resolve() if args.outdir else Path(__file__).resolve().parent

    if args.mode == "demo":
        build_demo_outputs(project_dir, peak)
        print(f"Demo outputs written to: {project_dir}")
        run_solve_mode(peak, args.lambda_values)
    elif args.mode == "solve":
        run_solve_mode(peak, args.lambda_values)
    elif args.mode == "curve":
        run_curve_mode(peak, args.lambda_values, args.curve_points, args.curve_max_factor, project_dir)


if __name__ == "__main__":
    main()
