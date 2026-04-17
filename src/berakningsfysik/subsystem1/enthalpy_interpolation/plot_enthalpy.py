from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from interpolation import build_enthalpy_dataframe, RedlichKisterModel


def save_dataset(df, csv_path: Path) -> None:
    """Sparar dataframe till CSV."""
    df.to_csv(csv_path, index=False)
    print(f"Sparade dataset till: {csv_path}")

def plot_enthalpy_d1_d2(x_data, y_data, x_grid, y_grid, d1, d2, order, save_path: Path) -> None:
    """Plottning av blandningsentalpin och koncentrationsderivator"""
    fig, (ax1, ax2, ax3) = plt.subplots(3, 1, figsize=(7, 9), sharex=True)

    # Entalpi med datapunkter
    ax1.plot(x_grid, y_grid, "k-", linewidth=2, label="RK")
    ax1.plot(x_data, y_data, "bo", markersize=4, label="Datapunkter")
    ax1.set_ylabel(r"$\Delta H_{\text{mix}}$ (eV/atom)")
    ax1.set_title(f"Redlich–Kister-anpassning (ordning {order})")
    ax1.grid(True, alpha=0.3)

    # Förstaderivata
    ax2.plot(x_grid, d1, "r-", linewidth=2)
    ax2.set_ylabel(r"$d\Delta H_{\text{mix}}/dx$ (eV/atom)")
    ax2.axhline(0, color="gray", linestyle=":", alpha=0.5)
    ax2.grid(True, alpha=0.3)

    # Andraderivata
    ax3.plot(x_grid, d2, "b--", linewidth=2)
    ax3.set_xlabel("Al-sammansättning $x$")
    ax3.set_ylabel(r"$d^2\Delta H_{\text{mix}}/dx^2$ (eV/atom)")
    ax3.axhline(0, color="gray", linestyle=":", alpha=0.5)
    ax3.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig(save_path, dpi=200)
    plt.close()
    print(f"Sparade figur till: {save_path}")


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Plotta interpolerad blandningsentalpi och dess derivator från QE .out-filer."
    )
    parser.add_argument(
        "data_dir",
        nargs="?",
        default="qe_outputs",
        help="Katalog med .out-filer. Standard: qe_outputs",
    )
    parser.add_argument(
        "--glob",
        default="*x=*.out",
        help="Glob-mönster för filnamn. Standard: *x=*.out",
    )
    parser.add_argument(
        "--order",
        type=int,
        default=3,
        help="Ordning på Redlich-Kister-polynomet. Standard: 3",
    )
    args = parser.parse_args()

    # Läs data och anpassa modell

    df = build_enthalpy_dataframe(args.data_dir, glob_pattern=args.glob)
    x = df["x"].to_numpy()
    y = df["H_mix_eV_per_atom"].to_numpy()
    model = RedlichKisterModel.fit(x=x, hmix=y, order=args.order)

    # Beräkna kurvor på grid
    x_grid = np.linspace(0.0, 1.0, 500)
    y_grid = model.hmix(x_grid)
    d1_grid = model.d1(x_grid)
    d2_grid = model.d2(x_grid)

    H_derivative = model.d1(x_grid)
    H_second_derivative = model.d2(x_grid)

    out_dir = Path(args.data_dir)
    csv_path = out_dir / "enthalpy_dataset.csv"
<<<<<<< HEAD
    fig_path = out_dir / "interpolated_hmix.png"

    df.to_csv(csv_path, index=False)

    plt.figure(figsize=(7, 4.5))
    plt.plot(x, y, "o", label="Datapunkter")
    plt.plot(x_grid, y_grid, "-", label=f"Redlich–Kister-fit (ordning {args.order})")
    plt.plot(x_grid, H_second_derivative, "-", label=f"Andraderivata av Redlich–Kister-fit (ordning {args.order})")
    plt.plot(x_grid, H_derivative, "-", label=f"Derivata av Redlich–Kister-fit (ordning {args.order})")
    plt.xlabel("Al-sammansättning x")
    plt.ylabel(r"Blandningsentalpi $\Delta H_{mix}$ (eV/atom)")
    plt.title("Interpolerad blandningsentalpi")
    plt.grid(True, alpha=0.3)
    plt.legend()
    plt.tight_layout()
    plt.savefig(fig_path, dpi=200)
    plt.close()
=======
    save_dataset(df, csv_path)
>>>>>>> 2691012e0bef1bb87d7ec54d2694a25967ce79c1

    # Skriv ut modellinfo
    print(f"Läste {len(df)} filer från: {out_dir.resolve()}")
    print("Redlich-Kister-koefficienter L_i [eV/atom]:")
    for i, c in enumerate(model.coeffs):
        print(f"L{i} = {c:.8f}")
    print(f"RMSE = {model.rmse:.6e} eV/atom")

    # Figurer
    fig_path = out_dir / "enthalpy_w_derivatives.png"
    plot_enthalpy_d1_d2(x, y, x_grid, y_grid, d1_grid, d2_grid, args.order, fig_path)

if __name__ == "__main__":
    main()