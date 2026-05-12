"""Denna fil sparar interpolerad deltaH, deltaS och deltaG samt
andraderivatan av deltaG (för alla temp)
i en csv-fil som enkelt kan läsas till en pandas dataframe med kommando
pd.read_csv("dataframe.csv")"""

import pandas as pd
import numpy as np
import argparse
from pathlib import Path

from enthalpy_interpolation.interpolation import RedlichKisterModel
from thermodynamics import entropy_per_atom
from spinodal_functions import entropy_second_derivative

DEFAULT_SYSTEM = "TiAlN"
DEFAULT_RESULTS_DIRNAME = "results"
DEFAULT_MODEL_DIRNAME = "rk_model"
DEFAULT_THERMODYNAMICS_DIRNAME = "thermodynamics"
DEFAULT_DATAFRAME_NAME = "dataframe.csv"

columns = ["x", "deltaH", "deltaS"]
temps = [i for i in range(0, 15001, 500)]
num_of_inter_points = 500
n = 2


def write_file():

    parser = argparse.ArgumentParser(description="Skapa termodynamisk dataframe.")
    parser.add_argument(
        "--system",
        default=DEFAULT_SYSTEM,
        help="Materialsystem. Standard: TiAlN.",
    )
    parser.add_argument(
        "--model-path",
        default=None,
        help=(
            "Sökväg till sparad RK-modell eller modellmapp. "
            "Standard: results/<system>/rk_model."
        ),
    )
    parser.add_argument(
        "--output",
        default=None,
        help=(
            "Sökväg till dataframe.csv. "
            "Standard: results/<system>/thermodynamics/dataframe.csv."
        ),
    )

    args = parser.parse_args()

    model_path = Path(args.model_path) if args.model_path else default_model_dir(args.system)
    output_path = Path(args.output) if args.output else default_dataframe_path(args.system)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    print(f"Räknar ut termodynamiska storheter.")
    print(f"Läser RK-modell från: {model_path}")
    print(f"Skriver data till: {output_path}")

    x_interpolation = np.linspace(0.0, 1.0, num_of_inter_points)
    model = find_model(model_path)
    H_interpolated = find_deltaH(model, x_interpolation)

    # skriver alla kolumner samt interpolationspunkter
    with open(output_path, "w", newline="") as file:

        file.write(columns[0])

        for column in columns[1:]:
            file.write(","+column)

        for i in range(len(temps)):
            file.write(",deltaG_T"+str(i))

        file.write("\n")

        for x in x_interpolation:
            file.write(str(x)+",\n")

    file.close()

    write_data(output_path, "deltaH", H_interpolated)
    write_data(output_path, "deltaS", find_deltaS(2, x_interpolation))
    

    for i in range(len(temps)):
        write_data(output_path, "deltaG_T" + str(i), find_deltaG(output_path, temps[i]))
        write_data(output_path, "d2deltaG_T" + str(i), find_d2G(model, x_interpolation, temps[i], n))
        

def write_data(csv_path, column, data_array):
    """Skriver datan i "data_array"
    under "column" i csv-filen.
    column = sträng från columns (t.ex "deltaH")
    data = vektor med data till kolumnen"""
    df = pd.read_csv(csv_path)
    df[column] = data_array
    df.to_csv(csv_path, index=False)

def find_deltaS(n, x_interpolated):
    """Hittar interpolerad deltaS_mix 
    n = atomer per metallplats"""
    deltaS = []
    for x in x_interpolated:
        deltaS.append(entropy_per_atom(x, n))
    return np.array(deltaS)

def find_deltaG(csv_path, T):
    df = pd.read_csv(csv_path)
    delta_S = df["deltaS"]
    delta_H = df["deltaH"]
    deltaG = delta_H - T*delta_S
    return deltaG

def load_saved_model(model_path):
    """
    Läser en sparad Redlich-Kister-modell.

    Modellen finns i subsystem_python/results/TiAlN/rk_model/

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
        data = np.load(model_path)
        coeffs = data["coeffs"]
        rmse = float(data["rmse"]) if "rmse" in data else np.nan
        return RedlichKisterModel(coeffs=coeffs, rmse=rmse)

    if model_path.suffix == ".npy":
        coeffs = np.load(model_path)
        return RedlichKisterModel(coeffs=coeffs, rmse=np.nan)

    raise ValueError(
        f"Okänd modellfiltyp: {model_path.suffix}. "
        "Använd .npz eller .npy."
    )


def find_model(model_path):
    return load_saved_model(model_path)

def find_deltaH(model, x_interpolation):
    """Hämtar interpolerade deltaH_mix
    x_interpolation = array med interpolationspunkter (använd np.linspace(...))"""
    
    return model.hmix(x_interpolation)

def find_d2G(model, x_interpolation, T, n):
    """returnerar vektor med 
    andraderivatan av deltaG"""
    d2deltaH = model.d2(x_interpolation)[1:-1]
    d2deltaS = []
    for x in x_interpolation[1:-1]:
        d2deltaS.append(entropy_second_derivative(x, n))

    d2deltaG = np.insert((d2deltaH - T*np.array(d2deltaS)), 0, np.nan)
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
