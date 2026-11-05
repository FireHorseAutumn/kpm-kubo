#!/usr/bin/env python3
"""R4 clean-metal figure: KPM sigma_1(omega) at L=128, quarter filling,
overlaid with the exact Lehmann spectrum of the separable open lattice
(binned) at the same e_f. The clean metal's DC content is the Drude
weight at omega = 0 (annotated, gated in test_r4); the stated resolution
Delta = 4 kernel half-widths is marked. The finite regular tail is the
open-boundary finite-L physics (~1/L), not noise. Reads
data/sigma_clean_metal_L128{,_exact}.csv, writes
data/sigma_clean_metal_L128.png."""

import pathlib

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

DATA = pathlib.Path(__file__).resolve().parent.parent / "data"

# the gated config (test_r4); Delta = 4 pi a / n_cheb
A = 4.04
N_CHEB = 128
DELTA = 4.0 * np.pi * A / N_CHEB


def load_csv(path):
    return np.genfromtxt(path, delimiter=",", skip_header=1)


def main():
    kpm = load_csv(DATA / "sigma_clean_metal_L128.csv")
    exact = load_csv(DATA / "sigma_clean_metal_L128_exact.csv")

    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    for ax, ylim in ((axes[0], None), (axes[1], 0.12)):
        ax.plot(kpm[:, 0], kpm[:, 1], lw=0.9,
                label="KPM (L=128, N_cheb=128, R=32, Jackson)")
        dwb = exact[1, 0] - exact[0, 0]
        ax.bar(exact[:, 0], exact[:, 1], width=dwb, alpha=0.45,
               label="exact Lehmann (separable spectrum)")
        ax.axvline(DELTA, color="k", ls=":", lw=1.0,
                   label=f"stated resolution $\\Delta$={DELTA:.3f}")
        ax.set_xlabel("$\\omega$")
        ax.set_xlim(0, 8)
        ax.legend()
    axes[0].set_ylabel("$\\sigma_1(\\omega)$")
    axes[0].set_title("clean 2D square-lattice metal, quarter filling, open BCs")
    axes[1].set_ylim(0, ylim)
    axes[1].set_title("regular part (zoom): finite-L weight ~1/L, not zero")
    fig.tight_layout()
    out = DATA / "sigma_clean_metal_L128.png"
    fig.savefig(out, dpi=150)
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
