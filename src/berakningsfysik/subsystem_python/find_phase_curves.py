"""Denna fil hittar kompositioner som gäller vid binodal- och
spinodal-kurvorna, och skriver dessa till curves.csv."""

from scipy.spatial import ConvexHull
from pathlib import Path
import argparse

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

from help_functions import find_parameters
from thermodynamics import T_string


DEFAULT_SYSTEM = "TiAlN"
DEFAULT_RESULTS_DIRNAME = "results"
DEFAULT_THERMODYNAMICS_DIRNAME = "thermodynamics"
DEFAULT_PHASE_CURVES_DIRNAME = "phase_curves"
DEFAULT_GIBBS_PLOTS_DIRNAME = "Gibbs_plots"
DEFAULT_DATAFRAME_NAME = "dataframe.csv"
DEFAULT_CURVES_NAME = "curves.csv"

columns = ["T", "xa", "xb", "spinodal_xa", "spinodal_xb"]


def main():
    parser = argparse.ArgumentParser(
        description="Hitta binodal- och spinodalkurvor."
    )

    parser.add_argument(
        "alloy_name",
        nargs="?",
        default=DEFAULT_SYSTEM,
        help="Legeringsnamn, t.ex. TiAlN. Parameterfil läses från alloy_parameters/<legering>.csv.",
    )

    parser.add_argument(
        "--input",
        default=None,
        help=(
            "Sökväg till dataframe.csv. "
            "Standard: results/<legering>/thermodynamics/dataframe.csv."
        ),
    )

    parser.add_argument(
        "--output",
        default=None,
        help=(
            "Sökväg till curves.csv. "
            "Standard: results/<legering>/phase_curves/curves.csv."
        ),
    )

    parser.add_argument(
        "--gibbs-plot-dir",
        default=None,
        help=(
            "Mapp för Gibbs/common-tangent-plottar. "
            "Standard: results/<legering>/phase_curves/Gibbs_plots."
        ),
    )

    args = parser.parse_args()

    # Läs temperaturer från alloy_parameters/<legering>.csv.
    _, temps, alloy_name, _ = find_parameters(args.alloy_name)

    input_path = (
        Path(args.input)
        if args.input
        else default_dataframe_path(alloy_name)
    )

    output_path = (
        Path(args.output)
        if args.output
        else default_curves_path(alloy_name)
    )

    gibbs_plot_dir = (
        Path(args.gibbs_plot_dir)
        if args.gibbs_plot_dir
        else default_gibbs_plot_dir(alloy_name)
    )

    output_path.parent.mkdir(parents=True, exist_ok=True)
    gibbs_plot_dir.mkdir(parents=True, exist_ok=True)

    print(f"Letar efter binodal- och spinodalkurvor för {alloy_name}.")
    print(f"Läser termodynamisk dataframe från: {input_path}")
    print(f"Sparar kurvor till: {output_path}")
    print(f"Sparar Gibbs-plottar till: {gibbs_plot_dir}")

    df = pd.read_csv(input_path)

    check_dataframe_has_temperature_columns(df, temps)

    with open(output_path, "w", newline="") as file:
        file.write(columns[0])

        for column in columns[1:]:
            file.write("," + column)

        file.write("\n")

        for T in temps:
            file.write(str(T) + ",\n")

    xa_list = []
    xb_list = []
    spinodal_xa_list = []
    spinodal_xb_list = []

    for index, T in enumerate(temps):
        xa, xb, spinodal_xa, spinodal_xb = find_comps_at_temp(
            T=T,
            df=df,
            index=index,
            plot_dir=gibbs_plot_dir,
        )

        xa_list.append(xa)
        xb_list.append(xb)
        spinodal_xa_list.append(spinodal_xa)
        spinodal_xb_list.append(spinodal_xb)

        print(
            f"T = {T:8}: "
            f"binodal = ({xa}, {xb}), "
            f"spinodal = ({spinodal_xa}, {spinodal_xb})"
        )

    df_curves = pd.read_csv(output_path)

    df_curves["xa"] = xa_list
    df_curves["xb"] = xb_list
    df_curves["spinodal_xa"] = spinodal_xa_list
    df_curves["spinodal_xb"] = spinodal_xb_list

    df_curves.to_csv(output_path, index=False)


def check_dataframe_has_temperature_columns(df: pd.DataFrame, temps: list) -> None:
    """
    Kontrollerar att dataframe.csv innehåller de deltaG- och d2deltaG-kolumner
    som motsvarar temperaturerna från parameterfilen.
    """
    missing_columns = []

    for i in range(len(temps)):
        for column in [f"deltaG_T{i}", f"d2deltaG_T{i}"]:
            if column not in df.columns:
                missing_columns.append(column)

    if missing_columns:
        raise ValueError(
            "dataframe.csv saknar kolumner som behövs för phase curve-beräkningen.\n"
            f"Saknade kolumner: {missing_columns}\n"
            "Kontrollera att build_dataframe.py kördes med samma parameterfil/temperaturer."
        )


def find_spinodal_points(x, d2G):
    """
    Hittar spinodalpunkter där d2G byter tecken.
    Använder linjär interpolation mellan gridpunkter.
    """
    spinodals = []

    x = np.asarray(x, dtype=float)
    d2G = np.asarray(d2G, dtype=float)

    for i in range(len(x) - 1):
        y0 = d2G[i]
        y1 = d2G[i + 1]

        if not np.isfinite(y0) or not np.isfinite(y1):
            continue

        if y0 == 0:
            spinodals.append(float(x[i]))

        elif y0 * y1 < 0:
            x0 = x[i]
            x1 = x[i + 1]

            # Linjär interpolation av nollstället:
            # y = 0 mellan (x0, y0) och (x1, y1).
            x_zero = x0 - y0 * (x1 - x0) / (y1 - y0)
            spinodals.append(float(x_zero))

    return spinodals


def find_lower_hull_binodal(x, deltaG, tol=1e-10):
    """
    Hittar binodalen från undre konvexa höljet.

    Returnerar xa, xb.
    Om ingen common tangent hittas returneras np.nan, np.nan.
    """
    x = np.asarray(x, dtype=float)
    deltaG = np.asarray(deltaG, dtype=float)

    finite_mask = np.isfinite(x) & np.isfinite(deltaG)
    x = x[finite_mask]
    deltaG = deltaG[finite_mask]

    points = np.array(list(zip(x, deltaG)))
    hull = ConvexHull(points)

    candidates = []

    for simplex in hull.simplices:
        i, j = sorted(simplex)

        # Intilliggande punkter är bara vanliga kurvsegment, inte en tie-line.
        if j - i <= 1:
            continue

        x_i = x[i]
        x_j = x[j]
        g_i = deltaG[i]
        g_j = deltaG[j]

        # Linje mellan två hullpunkter.
        line = g_i + (g_j - g_i) * (x - x_i) / (x_j - x_i)

        # Undre hull: alla punkter ska ligga ovanför eller på linjen.
        if np.all(deltaG >= line - tol):
            width = x_j - x_i
            candidates.append((width, x_i, x_j))

    if not candidates:
        return np.nan, np.nan

    # Om flera kandidater finns: ta den bredaste blandbarhets gapen.
    candidates.sort(reverse=True)
    _, xa, xb = candidates[0]

    return float(xa), float(xb)


def find_comps_at_temp(T, df, index, plot_dir):
    """
    Returnerar kompositioner för binodal och spinodal vid temperaturen T.

    Returnerar:
        xa, xb, spinodal_xa, spinodal_xb

    Alla är floats eller np.nan.
    """
    x_interpolated = df["x"].to_numpy(dtype=float)
    deltaG = df["deltaG_T" + str(index)].to_numpy(dtype=float)
    d2deltaG = df["d2deltaG_T" + str(index)].to_numpy(dtype=float)

    # Hitta spinodalpunkter från d2G = 0.
    spinodals = find_spinodal_points(x_interpolated, d2deltaG)

    # Hitta binodal från undre konvexa höljet.
    xa, xb = find_lower_hull_binodal(x_interpolated, deltaG)

    if np.isfinite(xa) and np.isfinite(xb):
        # Ta spinodalpunkter som ligger mellan binodalpunkterna.
        inner_spinodals = [
            s for s in spinodals
            if xa < s < xb
        ]

        if len(inner_spinodals) >= 2:
            spinodal_xa = inner_spinodals[0]
            spinodal_xb = inner_spinodals[-1]
        else:
            spinodal_xa = np.nan
            spinodal_xb = np.nan

    else:
        # Ingen binodal hittades.
        if len(spinodals) >= 2:
            xa = np.nan
            xb = np.nan
            spinodal_xa = spinodals[0]
            spinodal_xb = spinodals[-1]
        else:
            xa = np.nan
            xb = np.nan
            spinodal_xa = np.nan
            spinodal_xb = np.nan

    plot_gibbs_with_common_tangent(
        x_interpolated=x_interpolated,
        deltaG=deltaG,
        T=T,
        xa=xa,
        xb=xb,
        spinodal_xa=spinodal_xa,
        spinodal_xb=spinodal_xb,
        plot_dir=plot_dir,
    )

    return xa, xb, spinodal_xa, spinodal_xb


def plot_gibbs_with_common_tangent(
    x_interpolated,
    deltaG,
    T,
    xa,
    xb,
    spinodal_xa,
    spinodal_xb,
    plot_dir,
):
    """
    Plottar deltaG, eventuell gemensam tangent och eventuella spinodalpunkter.
    """
    plt.clf()
    plt.plot(x_interpolated, deltaG, "b-", label=r"$\Delta G_{mix}$")

    if np.isfinite(xa) and np.isfinite(xb):
        g_xa = np.interp(xa, x_interpolated, deltaG)
        g_xb = np.interp(xb, x_interpolated, deltaG)

        plt.plot([xa, xb], [g_xa, g_xb], "k--", label="gemensam tangent")
        plt.plot([xa, xb], [g_xa, g_xb], "ko")

    if np.isfinite(spinodal_xa):
        plt.axvline(spinodal_xa, color="r", linestyle=":", alpha=0.7)

    if np.isfinite(spinodal_xb):
        plt.axvline(spinodal_xb, color="r", linestyle=":", alpha=0.7)

    plt.xlabel("x")
    plt.ylabel(r"$\Delta G_{mix}$")
    plt.title(f"T = {T} K")
    plt.legend()
    plt.grid(True, alpha=0.3)

    plot_dir = Path(plot_dir)
    plot_dir.mkdir(parents=True, exist_ok=True)

    plt.savefig(plot_dir / f"hull_deltaG_T={T_string(T)}.png", dpi=200)


def project_root() -> Path:
    """
    Returnerar subsystem_python-mappen.
    """
    return Path(__file__).resolve().parent


def default_results_root() -> Path:
    return project_root() / DEFAULT_RESULTS_DIRNAME


def default_system_dir(system: str) -> Path:
    return default_results_root() / system


def default_dataframe_path(system: str) -> Path:
    return (
        default_system_dir(system)
        / DEFAULT_THERMODYNAMICS_DIRNAME
        / DEFAULT_DATAFRAME_NAME
    )


def default_phase_curves_dir(system: str) -> Path:
    return default_system_dir(system) / DEFAULT_PHASE_CURVES_DIRNAME


def default_curves_path(system: str) -> Path:
    return default_phase_curves_dir(system) / DEFAULT_CURVES_NAME


def default_gibbs_plot_dir(system: str) -> Path:
    return default_phase_curves_dir(system) / DEFAULT_GIBBS_PLOTS_DIRNAME


if __name__ == "__main__":
    main()