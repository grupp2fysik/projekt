#!/usr/bin/env python3
from __future__ import annotations

import argparse
import sys
from pathlib import Path

# Lägg till subsystem_python/ i Python-sökvägen.
# Denna fil ligger i subsystem_python/validation/scripts/,
# så parents[2] är subsystem_python/.
SUBSYSTEM_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(SUBSYSTEM_ROOT))

import numpy as np
import matplotlib.pyplot as plt

from enthalpy_interpolation.interpolation import (
    build_enthalpy_dataframe,
    RedlichKisterModel,
)

KB_EV_PER_K = 8.6173e-5
TEMPERATURES = [2000, 4000, 6000, 8000]


def main() -> None:
    parser = argparse.ArgumentParser(
        description=(
            "Plotta interpolerad blandningsentalpi och Gibbs fria blandningsenergi "
            "vid flera temperaturer i eV/formelenhet."
        )
    )

    parser.add_argument(
        "data_csv",
        help="CSV-fil med x och blandningsentalpi.",
    )

    parser.add_argument(
        "--model-path",
        required=True,
        help="Sökväg till sparad RK-modell, t.ex. rk_model/rk_model.npz.",
    )

    parser.add_argument(
        "--hmix-column",
        default="hmix_ev_per_formula_unit",
        help="Kolumn för blandningsentalpi i inputfilen.",
    )

    parser.add_argument(
        "--atoms-per-formula-unit",
        type=float,
        default=2.0,
        help="Antal atomer per formelenhet. För TiAlN är detta 2.",
    )

    parser.add_argument(
        "--num-points",
        type=int,
        default=600,
        help="Antal punkter i interpolerad kurva.",
    )

    parser.add_argument(
        "--output",
        required=True,
        help="Sökväg till figuren som ska sparas.",
    )

    args = parser.parse_args()

    data_csv = Path(args.data_csv)
    model_path = Path(args.model_path)
    output_path = Path(args.output)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    df = build_enthalpy_dataframe(
        data_csv,
        hmix_column=args.hmix_column,
        atoms_per_formula_unit=args.atoms_per_formula_unit,
    )

    model = load_rk_model(model_path)

    x_data = df["x"].to_numpy(dtype=float)
    h_data_fu = df["H_mix_eV_per_formula_unit"].to_numpy(dtype=float)

    x_grid = np.linspace(0.0, 1.0, args.num_points)

    h_grid_fu = model.hmix(x_grid)

    fig, ax = plt.subplots(figsize=(7.5, 5.2))

    # Tydliga pastellfärger.
    color_hmix = "#7B68EE"      # pastell-lila
    color_points = "#222222"    # nästan svart, tydliga datapunkter

    gibbs_colors = {
        2000: "#66C2A5",  # pastell-grön
        4000: "#8DA0CB",  # pastell-blå/lavendel
        6000: "#FC8D62",  # pastell-orange
        8000: "#E78AC3",  # pastell-rosa
    }

    ax.plot(
        x_grid,
        h_grid_fu,
        color=color_hmix,
        linewidth=2.8,
        label=r"$\Delta H_{\mathrm{mix}}$",
    )

    ax.plot(
        x_data,
        h_data_fu,
        "o",
        color=color_points,
        markerfacecolor="white",
        markeredgewidth=1.4,
        markersize=5.5,
        label="Datapunkter",
    )

    for T in TEMPERATURES:
        g_grid_fu = h_grid_fu - T * entropy_formula_unit(x_grid)
        ax.plot(
            x_grid,
            g_grid_fu,
            color=gibbs_colors[T],
            linewidth=2.4,
            label=rf"$\Delta G_{{\mathrm{{mix}}}}$, {T} K",
        )

    ax.axhline(0.0, color="#555555", linestyle=":", linewidth=1.2)

    ax.set_xlabel(r"$x$ i $\mathrm{Ti}_{1-x}\mathrm{Al}_x\mathrm{N}$")
    ax.set_ylabel(r"Energi (eV/formelenhet)")
    ax.set_title(
        rf"Interpolerad $\Delta H_{{\mathrm{{mix}}}}$ och "
        rf"$\Delta G_{{\mathrm{{mix}}}}$"
    )

    ax.set_xlim(0.0, 1.0)
    ax.grid(True, alpha=0.3)
    ax.legend(frameon=True)

    fig.tight_layout()
    fig.savefig(output_path, dpi=250)
    plt.close(fig)

    print(f"Sparade sammansatt Hmix/Gibbs-figur till: {output_path}")


def entropy_formula_unit(x: np.ndarray) -> np.ndarray:
    """
    Ideal konfigurationell blandningsentropi per formelenhet.

    För Ti_{1-x}Al_xN sker blandningen på en metallplats per formelenhet.
    Därför är:

        S_mix = -kB [x ln x + (1-x) ln(1-x)]

    Enhet: eV/(formelenhet K)
    """
    x = np.asarray(x, dtype=float)
    S = np.zeros_like(x)

    mask = (x > 0.0) & (x < 1.0)

    xm = x[mask]
    S[mask] = -KB_EV_PER_K * (
        xm * np.log(xm) + (1.0 - xm) * np.log(1.0 - xm)
    )

    return S


def load_rk_model(model_path: Path) -> RedlichKisterModel:
    """
    Läser RK-modell från .npz eller .npy.
    """
    if model_path.is_dir():
        model_path = model_path / "rk_model.npz"

    if not model_path.exists():
        raise FileNotFoundError(f"Hittar inte RK-modellen: {model_path}")

    if model_path.suffix == ".npz":
        with np.load(model_path, allow_pickle=False) as data:
            coeffs = np.asarray(data["coeffs"], dtype=float)
            rmse = float(data["rmse"]) if "rmse" in data else np.nan

        return RedlichKisterModel(coeffs=coeffs, rmse=rmse)

    if model_path.suffix == ".npy":
        coeffs = np.asarray(np.load(model_path), dtype=float)
        return RedlichKisterModel(coeffs=coeffs, rmse=np.nan)

    raise ValueError(f"Okänd modellfiltyp: {model_path.suffix}")


if __name__ == "__main__":
    main()