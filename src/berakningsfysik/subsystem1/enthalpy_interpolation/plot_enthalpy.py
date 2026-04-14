from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from interpolation import build_enthalpy_dataframe, RedlichKisterModel


def main() -> None:
    parser = argparse.ArgumentParser(description="Plotta interpolerad blandningsentalpi från QE .out-filer.")
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

    df = build_enthalpy_dataframe(args.data_dir, glob_pattern=args.glob)

    x = df["x"].to_numpy()
    y = df["H_mix_eV_per_atom"].to_numpy()

    model = RedlichKisterModel.fit(x=x, hmix=y, order=args.order)

    x_grid = np.linspace(0.0, 1.0, 500)
    y_grid = model.hmix(x_grid)

    H_derivative = model.d1(x_grid)
    H_second_derivative = model.d2(x_grid)

    out_dir = Path(args.data_dir)
    csv_path = out_dir / "enthalpy_dataset.csv"
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

    print(f"Läste {len(df)} filer från: {out_dir.resolve()}")
    print(f"Sparade CSV till: {csv_path}")
    print(f"Sparade figur till: {fig_path}")
    print("Redlich-Kister-koefficienter [eV/atom]:")
    for i, c in enumerate(model.coeffs):
        print(f"L{i} = {c:.8f}")
    print(f"RMSE = {model.rmse:.6e} eV/atom")


if __name__ == "__main__":
    main()
