from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from interpolation import build_enthalpy_dataframe, RedlichKisterModel


DEFAULT_SYSTEM = "TiAlN"
DEFAULT_MODEL_DIRNAME = "rk_model"
DEFAULT_PLOT_DIRNAME = "interpolation_plot"
DEFAULT_MODEL_NPZ_NAME = "rk_model.npz"
DEFAULT_COEFFS_NPY_NAME = "rk_coeffs.npy"
DEFAULT_FIGURE_NAME = "enthalpy_w_derivatives.png"
DEFAULT_DATASET_NAME = "enthalpy_dataset.csv"


def project_root() -> Path:
    """
    Returnerar katalogen där denna fil ligger.
    Exempel:
        subsystem_python/enthalpy_interpolation/
    """
    return Path(__file__).resolve().parent


def subsystem_root() -> Path:
    """
    Returnerar huvudmappen subsystem_python/.
    """
    return project_root().parent


def default_results_root() -> Path:
    """
    Returnerar resultatmappen:
        subsystem_python/results/
    """
    return subsystem_root() / "results"


def default_system_dir(system: str) -> Path:
    """
    Returnerar systemets resultatmapp, t.ex.
        subsystem_python/results/TiAlN/
    """
    return default_results_root() / system


def default_rk_model_dir(system: str) -> Path:
    """
    Returnerar standardmappen för RK-modeller.
    """
    return default_system_dir(system) / DEFAULT_MODEL_DIRNAME


def default_interpolation_plot_dir(system: str) -> Path:
    return default_system_dir(system) / DEFAULT_PLOT_DIRNAME


def default_model_path(system: str) -> Path:
    """
    Standardmodellfilen.
    """
    return default_rk_model_dir(system) / DEFAULT_MODEL_NPZ_NAME


def _npz_scalar(npz_file, key: str, default=None):
    """
    Hämtar ett skalärt värde från en npz-fil.
    """
    if key not in npz_file:
        return default

    value = npz_file[key]

    if np.ndim(value) == 0:
        return value.item()

    return value


def load_rk_model(model_path: str | Path) -> tuple[RedlichKisterModel, dict]:
    """
    Läser in en sparad Redlich-Kister-modell.

    Föredrar ny .npz-fil:
        rk_model.npz

    men kan även läsa gammal .npy-fil:
        rk_coeffs.npy
    """
    model_path = Path(model_path)

    if not model_path.exists():
        raise FileNotFoundError(
            f"Hittar inte modellfilen {model_path}.\n"
            "Kör först interpolationprogrammet med --save_model."
        )

    metadata = {}

    if model_path.suffix == ".npz":
        with np.load(model_path, allow_pickle=False) as data:
            coeffs = np.asarray(data["coeffs"], dtype=float)
            rmse = float(_npz_scalar(data, "rmse", np.nan))

            metadata = {
                "rmse": rmse,
                "order": int(_npz_scalar(data, "order", len(coeffs) - 1)),
                "energy_basis": _npz_scalar(data, "energy_basis", None),
                "hmix_fit_col": _npz_scalar(data, "hmix_fit_col", None),
                "hmix_column": _npz_scalar(data, "hmix_column", None),
                "atoms_per_formula_unit": _npz_scalar(data, "atoms_per_formula_unit", None),
                "data_source": _npz_scalar(data, "data_source", None),
            }

        return RedlichKisterModel(coeffs=coeffs, rmse=rmse), metadata

    if model_path.suffix == ".npy":
        coeffs = np.asarray(np.load(model_path), dtype=float)
        metadata = {
            "rmse": np.nan,
            "order": len(coeffs) - 1,
            "energy_basis": None,
            "hmix_fit_col": None,
            "hmix_column": None,
            "atoms_per_formula_unit": None,
            "data_source": None,
        }
        return RedlichKisterModel(coeffs=coeffs, rmse=np.nan), metadata

    raise ValueError(
        f"Okänd modellfilstyp: {model_path.suffix}. "
        "Använd .npz eller .npy."
    )


def save_dataset(df, csv_path: Path) -> None:
    """Sparar dataframe till CSV."""
    df.to_csv(csv_path, index=False)
    print(f"Sparade dataset till: {csv_path}")


def plot_enthalpy_d1_d2(
    x_data: np.ndarray,
    y_data: np.ndarray,
    x_grid: np.ndarray,
    y_grid: np.ndarray,
    d1: np.ndarray,
    d2: np.ndarray,
    order: int,
    unit_label: str,
    save_path: Path,
) -> None:
    """Plottar blandningsentalpin och koncentrationsderivatorna."""
    fig, (ax1, ax2, ax3) = plt.subplots(3, 1, figsize=(7, 9), sharex=True)

    ax1.plot(x_grid, y_grid, "k-", linewidth=2, label="RK-modell")
    ax1.plot(x_data, y_data, "bo", markersize=4, label="Datapunkter")
    ax1.set_ylabel(rf"$\Delta H_{{\mathrm{{mix}}}}$ ({unit_label})")
    ax1.set_title(f"Redlich-Kister-modell från sparade koefficienter, ordning {order}")
    ax1.grid(True, alpha=0.3)
    ax1.legend()

    ax2.plot(x_grid, d1, "r-", linewidth=2)
    ax2.set_ylabel(rf"$d\Delta H_{{\mathrm{{mix}}}}/dx$ ({unit_label})")
    ax2.axhline(0, color="gray", linestyle=":", alpha=0.5)
    ax2.grid(True, alpha=0.3)

    ax3.plot(x_grid, d2, "b--", linewidth=2)
    ax3.set_xlabel(r"Al-sammansättning $x$")
    ax3.set_ylabel(rf"$d^2\Delta H_{{\mathrm{{mix}}}}/dx^2$ ({unit_label})")
    ax3.axhline(0, color="gray", linestyle=":", alpha=0.5)
    ax3.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig(save_path, dpi=200)
    plt.close(fig)
    print(f"Sparade figur till: {save_path}")


def main() -> None:
    parser = argparse.ArgumentParser(
        description=(
            "Plotta blandningsentalpi och derivator från redan sparade "
            "Redlich-Kister-koefficienter. Input kan vara CSV-fil eller QE-katalog."
        )
    )

    parser.add_argument(
        "data_source",
        nargs="?",
        default="qe_outputs/qe_outputs2",
        help="CSV-fil med x och H_mix, eller katalog med QE .out-filer.",
    )

    parser.add_argument(
        "--model-path",
        default=None,
        help=(
            "Sökväg till sparad RK-modell. "
            "Standard: results/<system>/rk_model/rk_model.npz."
            "t.ex. results/TiAlN/rk_model/rk_model.npz."
        ),
    )

    parser.add_argument(
        "--glob",
        default="*x=*.out",
        help="Glob-mönster om data_source är en QE-katalog. Standard: *x=*.out",
    )

    parser.add_argument(
        "--hmix-column",
        default="hmix_ev_per_atom",
        help=(
            "H_mix-kolumn om data_source är CSV. Exempel: hmix_ev_per_atom, "
            "hmix_mev_per_atom, hmix_ry_per_atom eller hmix_ev_per_formula_unit."
        ),
    )

    parser.add_argument(
        "--atoms-per-formula-unit",
        type=float,
        default=2.0,
        help="Antal atomer per formelenhet. För Ti_{1-x}Al_xN är detta 2. Standard: 2.",
    )

    parser.add_argument(
        "--energy-basis",
        choices=["atom", "formula_unit"],
        default="atom",
        help=(
            "Vilken energibas som ska plottas. Måste matcha modellen du sparade. "
            "'atom' ger eV/atom. 'formula_unit' ger eV/formelenhet. Standard: atom."
        ),
    )

    parser.add_argument(
        "--num-points",
        type=int,
        default=500,
        help="Antal punkter i den släta modellkurvan. Standard: 500.",
    )

    parser.add_argument(
        "--output",
        default=None,
        help=(
            f"Sökväg till figuren. Standard: {DEFAULT_FIGURE_NAME} "
            "i samma katalog som data_source."
        ),
    )

    parser.add_argument(
        "--no-save-dataset",
        action="store_true",
        help="Spara inte det normaliserade entalpidatasetet till enthalpy_dataset.csv.",
    )

    parser.add_argument(
        "--system",
        default=DEFAULT_SYSTEM,
        help="Materialsystem. Styr standardmappar. Standard: TiAlN.",
    )

    parser.add_argument(
        "--plot-dir",
        default=None,
        help=(
            "Mapp där figuren ska sparas. "
            "Standard: results/<system>/interpolation_plot, "
            "t.ex. results/TiAlN/interpolation_plot."
        ),
    )

    args = parser.parse_args()

    data_source = Path(args.data_source)

    plot_dir = Path(args.plot_dir) if args.plot_dir else default_interpolation_plot_dir(args.system)
    plot_dir.mkdir(parents=True, exist_ok=True)

    model_path = Path(args.model_path) if args.model_path else default_model_path(args.system)
    fig_path = Path(args.output) if args.output else plot_dir / DEFAULT_FIGURE_NAME
    csv_path = plot_dir / DEFAULT_DATASET_NAME

    # Läs datapunkterna.
    # Om data_source är en katalog parsas QE-filerna och H_mix byggs som tidigare.
    # Om data_source är en CSV används H_mix direkt från vald kolumn.
    df = build_enthalpy_dataframe(
        data_source,
        glob_pattern=args.glob,
        hmix_column=args.hmix_column,
        atoms_per_formula_unit=args.atoms_per_formula_unit,
    )

    hmix_column = {
        "atom": "H_mix_eV_per_atom",
        "formula_unit": "H_mix_eV_per_formula_unit",
    }[args.energy_basis]

    unit_label = {
        "atom": "eV/atom",
        "formula_unit": "eV/formelenhet",
    }[args.energy_basis]

    x_data = df["x"].to_numpy(dtype=float)
    y_data = df[hmix_column].to_numpy(dtype=float)

    # Läs den anpassade modellen.
    model, metadata = load_rk_model(model_path)

    x_grid = np.linspace(0.0, 1.0, args.num_points)
    y_grid = model.hmix(x_grid)
    d1_grid = model.d1(x_grid)
    d2_grid = model.d2(x_grid)

    if not args.no_save_dataset:
        save_dataset(df, csv_path)

    print(f"Läste {len(df)} datapunkter från: {data_source.resolve()}")
    print(f"Läste RK-modell från: {model_path.resolve()}")
    print(f"Plottar energibas: {hmix_column} ({unit_label})")

    if metadata.get("energy_basis") is not None:
        print(f"Modellen sparades med energy_basis: {metadata['energy_basis']}")
    else:
        print("OBS: modellfilen innehåller ingen metadata om energibas.")

    print(f"Redlich-Kister-ordning: {model.order}")
    print(f"RMSE = {model.rmse:.6e} {unit_label}")

    print(f"Redlich-Kister-koefficienter L_i [{unit_label}]:")
    for i, c in enumerate(model.coeffs):
        print(f"L{i} = {c:.8f}")

    plot_enthalpy_d1_d2(
        x_data=x_data,
        y_data=y_data,
        x_grid=x_grid,
        y_grid=y_grid,
        d1=d1_grid,
        d2=d2_grid,
        order=model.order,
        unit_label=unit_label,
        save_path=fig_path,
    )


if __name__ == "__main__":
    main()
