"""Denna fil sparar interpolerad deltaH, deltaS och deltaG samt
andraderivatan av deltaG för alla temperaturer.

Data sparas till:
    results/<legering>/thermodynamics/dataframe.csv

Filen använder en redan sparad Redlich-Kister-modell från:
    results/<legering>/rk_model/rk_model.npz

Alltså: ingen ny interpolation görs här.
"""

import sys
import pandas as pd
import numpy as np
import argparse
from pathlib import Path

from enthalpy_interpolation.interpolation import RedlichKisterModel
from thermodynamics import entropy_per_atom, check_if_valid_T, find_temp_limits, check_if_valid_n
from spinodal_functions import entropy_second_derivative
from help_functions import find_parameters
from parameters import *


DEFAULT_SYSTEM = "TiAlN"
DEFAULT_RESULTS_DIRNAME = "results"
DEFAULT_MODEL_DIRNAME = "rk_model"
DEFAULT_THERMODYNAMICS_DIRNAME = "thermodynamics"
DEFAULT_DATAFRAME_NAME = "dataframe.csv"

columns = ["x", "deltaH", "deltaS"]


# Defaultvärden. Dessa skrivs över när vi läser alloy_parameters/<legering>.csv.
n = 2
temps = [i for i in range(0, 15001, 500)]
alloy_name = DEFAULT_SYSTEM


def _alloy_name_from_argv() -> str:
    """
    För kompatibilitet med andra filer som gör:
        from build_dataframe import temps

    Om scriptet körs som:
        python3 build_dataframe.py TiAlN

    eller om en annan fil körs som:
        python3 find_phase_curves.py TiAlN

    och importerar build_dataframe, försöker vi använda sys.argv[1]
    som legeringsnamn.
    """
    if len(sys.argv) > 1 and not sys.argv[1].startswith("-"):
        return sys.argv[1]

    return DEFAULT_SYSTEM


def _load_default_parameters_from_argv() -> None:
    """
    Försöker sätta globala n, temps och alloy_name från parameterfilen.

    Detta gör att gamla imports som
        from build_dataframe import temps
    fortfarande fungerar hyfsat med interfacet.
    """
    global n, temps, alloy_name

    candidate_alloy = _alloy_name_from_argv()

    try:
        n_read, temps_read, alloy_read, _ = find_parameters(candidate_alloy)
        n = n_read
        temps = temps_read
        alloy_name = alloy_read
    except Exception:
        # Vid t.ex. tester eller import utan parameterfil använder vi defaults.
        n = 2
        temps = [i for i in range(0, 15001, 500)]
        alloy_name = DEFAULT_SYSTEM


_load_default_parameters_from_argv()


def write_file():
    parser = argparse.ArgumentParser(description="Skapa termodynamisk dataframe.")
    parser.add_argument(
        "alloy_name",
        nargs="?",
        default=DEFAULT_SYSTEM,
        help="Legeringsnamn, t.ex. TiAlN. Parameterfil läses från alloy_parameters/<legering>.csv.",
    )
    parser.add_argument(
        "--model-path",
        default=None,
        help=(
            "Sökväg till sparad RK-modell eller modellmapp. "
            "Standard: results/<legering>/rk_model."
        ),
    )
    parser.add_argument(
        "--output",
        default=None,
        help=(
            "Sökväg till dataframe.csv. "
            "Standard: results/<legering>/thermodynamics/dataframe.csv."
        ),
    )

    args = parser.parse_args()

    n_local, temps_local, alloy_name_local, _ = find_parameters(args.alloy_name)

    model_path = (
        Path(args.model_path)
        if args.model_path
        else default_model_dir(alloy_name_local)
    )

    output_path = (
        Path(args.output)
        if args.output
        else default_dataframe_path(alloy_name_local)
    )

    output_path.parent.mkdir(parents=True, exist_ok=True)

    print(f"Räknar ut termodynamiska storheter för {alloy_name_local}.")
    print(f"Läser RK-modell från: {model_path}")
    print(f"Skriver data till: {output_path}")

    x_interpolation = np.linspace(0.0, 1.0, num_of_inter_points)
    model = find_model(model_path)
    H_interpolated = find_deltaH(model, x_interpolation)

    with open(output_path, "w", newline="") as file:
        file.write(columns[0])

        for column in columns[1:]:
            file.write("," + column)

        for i in range(len(temps_local)):
            file.write(",deltaG_T" + str(i))

        file.write("\n")

        for x in x_interpolation:
            file.write(str(x) + ",\n")

    write_data(output_path, "deltaH", H_interpolated)
    write_data(output_path, "deltaS", find_deltaS(n_local, x_interpolation))

    for i, T in enumerate(temps_local):
        write_data(output_path, "deltaG_T" + str(i), find_deltaG(output_path, T))
        write_data(
            output_path,
            "d2deltaG_T" + str(i),
            find_d2G(model, x_interpolation, T, n_local),
        )


def write_data(csv_path, column, data_array):
    """Skriver datan i data_array under column i csv-filen."""
    df = pd.read_csv(csv_path)
    df[column] = data_array
    df.to_csv(csv_path, index=False)


def find_deltaS(n, x_interpolated):
    """Hittar interpolerad deltaS_mix.

    n = atomer per metallplats / normaliseringsfaktor.
    """
    deltaS = []
    for x in x_interpolated:
        deltaS.append(entropy_per_atom(x, n))
    return np.array(deltaS)


def find_deltaG(csv_path, T):
    df = pd.read_csv(csv_path)
    delta_S = df["deltaS"]
    delta_H = df["deltaH"]
    deltaG = delta_H - T * delta_S
    return deltaG


def load_saved_model(model_path):
    """
    Läser en sparad Redlich-Kister-modell.

    Stödjer både:
      - rk_model.npz   där coeffs och rmse finns sparade
      - rk_coeffs.npy  där bara coeffs finns sparade
    """
    model_path = Path(model_path)

    if model_path.is_dir():
        npz_path = model_path / "rk_model.npz"
        npy_path = model_path / "rk_coeffs.npy"

        if npz_path.exists():
            model_path = npz_path
        elif npy_path.exists():
            model_path = npy_path
        else:
            raise FileNotFoundError(
                f"Hittar varken rk_model.npz eller rk_coeffs.npy i {model_path}"
            )

    if not model_path.exists():
        raise FileNotFoundError(f"Hittar inte modellfilen: {model_path}")

    if model_path.suffix == ".npz":
        with np.load(model_path, allow_pickle=False) as data:
            coeffs = np.asarray(data["coeffs"], dtype=float)
            rmse = float(data["rmse"]) if "rmse" in data else np.nan

        return RedlichKisterModel(coeffs=coeffs, rmse=rmse)

    if model_path.suffix == ".npy":
        coeffs = np.asarray(np.load(model_path), dtype=float)
        return RedlichKisterModel(coeffs=coeffs, rmse=np.nan)

    raise ValueError(
        f"Okänd modellfiltyp: {model_path.suffix}. "
        "Använd .npz eller .npy."
    )


def find_model(model_path):
    return load_saved_model(model_path)


def find_deltaH(model, x_interpolation):
    """Hämtar interpolerade deltaH_mix."""
    return model.hmix(x_interpolation)


def find_d2G(model, x_interpolation, T, n):
    """Returnerar vektor med andraderivatan av deltaG."""
    d2deltaH = model.d2(x_interpolation)[1:-1]

    d2deltaS = []
    for x in x_interpolation[1:-1]:
        d2deltaS.append(entropy_second_derivative(x, n))

    d2deltaG = np.insert((d2deltaH - T * np.array(d2deltaS)), 0, np.nan)
    return np.append(d2deltaG, [np.nan])


def project_root() -> Path:
    """
    Returnerar subsystem_python-mappen.
    build_dataframe.py ligger direkt i subsystem_python.
    """
    return Path(__file__).resolve().parent


def default_results_root() -> Path:
    """
    Returnerar subsystem_python/results.
    """
    return project_root() / DEFAULT_RESULTS_DIRNAME


def default_system_dir(system: str) -> Path:
    """
    Returnerar t.ex. subsystem_python/results/TiAlN.
    """
    return default_results_root() / system


def default_model_dir(system: str) -> Path:
    """
    Returnerar t.ex. subsystem_python/results/TiAlN/rk_model.
    """
    return default_system_dir(system) / DEFAULT_MODEL_DIRNAME


def default_dataframe_path(system: str) -> Path:
    """
    Returnerar t.ex. subsystem_python/results/TiAlN/thermodynamics/dataframe.csv.
    """
    return (
        default_system_dir(system)
        / DEFAULT_THERMODYNAMICS_DIRNAME
        / DEFAULT_DATAFRAME_NAME
    )


if __name__ == "__main__":
    write_file()