"""Generate Figure 6/7-style Weibull sensitivity plots.

Figure 6 confirms m = 1.8, 2.0, 2.2, 2.4, 2.6, but the supplied crop does
not print the fixed E, F0, or lambda. Figure 7's caption says F0 while its
legend says crack spacing (1.9 to 1.1 cm), so no exact F0 mapping is assumed.
The generated curves are therefore explicitly illustrative rather than
digitized traces of the published figures.
"""
from pathlib import Path

import matplotlib.pyplot as plt

from rock_weibull_model import PeakState, WeibullParameters, stress_strain_curve


def plot_family(peak, lambda_value, parameter_sets, out_stem, title):
    plt.figure(figsize=(6, 5))
    for label, params in parameter_sets:
        curve = stress_strain_curve(peak, lambda_value, params, n=500, F_max_factor=2.2)
        plt.plot(
            [epsilon for epsilon, _ in curve],
            [sigma for _, sigma in curve],
            label=label,
        )
    plt.xlabel("Axial strain")
    plt.ylabel("Axial stress (MPa)")
    plt.title(title)
    plt.xlim(left=0.0)
    plt.ylim(bottom=0.0)
    plt.grid(True, linestyle="--", alpha=0.35)
    plt.legend()
    plt.tight_layout()
    for suffix in (".png", ".svg"):
        plt.savefig(out_stem.with_suffix(suffix), dpi=200)
    plt.close()


def main():
    out_dir = Path(__file__).resolve().parent / "figures"
    out_dir.mkdir(parents=True, exist_ok=True)
    peak = PeakState(E=8500.0, nu=0.25, sigma3=0.0, sigma1c=30.0, epsilon1c=0.005, phi_deg=35.0)
    lambda_value = 0.90

    m_sets = [
        (f"m={m_value:.1f}", WeibullParameters(m=m_value, F0=25.0))
        for m_value in (1.8, 2.0, 2.2, 2.4, 2.6)
    ]
    plot_family(
        peak,
        lambda_value,
        m_sets,
        out_dir / "figure6_m_sensitivity",
        "Figure 6-style sensitivity to m (illustrative)",
    )

    f0_sets = [
        (f"F0={f0_value:.1f} MPa", WeibullParameters(m=2.2, F0=f0_value))
        for f0_value in (21.0, 23.0, 25.0, 27.0, 29.0)
    ]
    plot_family(
        peak,
        lambda_value,
        f0_sets,
        out_dir / "figure7_F0_sensitivity",
        "Figure 7-style sensitivity to F0 (illustrative)",
    )
    print(f"Wrote sensitivity figures to {out_dir}")


if __name__ == "__main__":
    main()
