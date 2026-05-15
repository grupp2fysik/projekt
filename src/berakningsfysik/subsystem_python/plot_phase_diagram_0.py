"""Denna fil använder datan i curves.csv för att plotta fasdiagrammet."""

from pathlib import Path
import argparse

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

from help_functions import find_parameters
from find_phase_curves import columns
from parameters import *


DEFAULT_SYSTEM = "TiAlN"
DEFAULT_RESULTS_DIRNAME = "results"
DEFAULT_MODEL_DIRNAME = "rk_model"
DEFAULT_THERMODYNAMICS_DIRNAME = "thermodynamics"
DEFAULT_DATAFRAME_NAME = "dataframe.csv"

def main():
    parser = argparse.ArgumentParser(description="Plotta fasdiagram.")
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
            "Sökväg till curves.csv. "
            "Standard: results/<legering>/phase_curves/curves.csv."
        ),
    )
    parser.add_argument(
        "--output",
        default=None,
        help=(
            "Sökväg till fasdiagrammet. "
            "Standard: results/<legering>/phase_curves/phase_diagram/phase_diagram.png."
        ),
    )

    args = parser.parse_args()

    #_, temps, alloy_name, _ = find_parameters(args.alloy_name)

    curves_path = (
        Path(args.input)
        if args.input
        else default_curves_path(alloy_name)
    )

    output_path = (
        Path(args.output)
        if args.output
        else default_phase_diagram_path(alloy_name)
    )

    output_path.parent.mkdir(parents=True, exist_ok=True)

    print(f"Plottar fasdiagram för {alloy_name}.")
    print(f"Läser kurvor från: {curves_path}")
    print(f"Sparar fasdiagram till: {output_path}")
    print(f"Sparar fasdiagram till: {plots_dirname}")

    df = pd.read_csv(curves_path)

    fig = plt.figure()

    plot_curve(df, temps, "binodal", columns[1], columns[2])
    plot_curve(df, temps, "spinodal", columns[3], columns[4])

    ticks = list(plt.yticks()[0])
    ticks += [temps[0], temps[-1]]
    ticks = sorted(set(ticks))
    plt.yticks(ticks)

    plt.xlim(0, 1)
    plt.ylim(temps[0], temps[-1])
    plt.xlabel("x")
    plt.ylabel("T [K]")
    plt.title(f"Fasdiagram {alloy_name}")
    fig.legend()

    plt.tight_layout()
    plt.savefig(output_path, dpi=200)
    plt.savefig(plots_dirname + f"/{alloy_name}_phase_diagram")
    plt.close(fig)


def plot_curve(df, temps, curve_type, column1, column2):
    """Sätter ihop kolumnerna column1 och column2 och plottar mot temperaturen.

    Exempel:
        column1 = xa
        column2 = xb

    curve_type är "spinodal" eller "binodal".
    """

    line_styles = {"spinodal": "r--", "binodal": "b-"}

    main_list1, max_num_elements1 = return_comps_main_list(df, column1)
    main_list2, max_num_elements2 = return_comps_main_list(df, column2)

    main_list2.reverse()
    max_num_elements = max(max_num_elements1, max_num_elements2)

    for j in range(0, max_num_elements):
        x_values = []
        y_values = []

        for index, comps_list in enumerate(main_list1):
            if len(comps_list) > j:
                if not np.isnan(comps_list[j]):
                    x_values.append(comps_list[j])
                    y_values.append(temps[index])

        for index, comps_list in enumerate(main_list2):
            if len(comps_list) > j:
                if not np.isnan(comps_list[j]):
                    x_values.append(comps_list[j])
                    y_values.append(temps[-index - 1])

        plt.plot(x_values, y_values, line_styles[curve_type], label=curve_type)


def turn_string_to_list(value):
    """
    Gör om ett värde från curves.csv till en lista med floats.

    Stödjer både gamla formatet:
        "[0.123, 0.456]"

    och nya formatet:
        0.123

    samt nan.
    """

    # Nya curves.csv: vanliga numeriska värden.
    if isinstance(value, (int, float, np.integer, np.floating)):
        if np.isnan(value):
            return [np.nan]
        return [float(value)]

    # Pandas-NaN eller None.
    if pd.isna(value):
        return [np.nan]

    # Gamla curves.csv: listor sparade som strängar.
    value = str(value).strip()

    if value == "" or value.lower() == "nan":
        return [np.nan]

    comps_string = value.strip("[]")

    if comps_string == "":
        return [np.nan]

    comps_list = comps_string.split(",")

    new_list = []
    for element in comps_list:
        element = element.strip()

        if element.lower() == "nan":
            new_list.append(np.nan)
        else:
            new_list.append(float(element))

    return new_list


def return_comps_main_list(df, column):
    """Tar en kolumn från curves.csv och omvandlar den till en nästlad lista."""

    comps_main_list = []
    max_num_elements = 1

    for index in df.index:
        comps_list = turn_string_to_list(df[column][index])
        comps_main_list.append(comps_list)

        if len(comps_list) > max_num_elements:
            max_num_elements = len(comps_list)

    return comps_main_list, max_num_elements


def project_root() -> Path:
    """
    Returnerar subsystem_python-mappen.
    plot_phase_diagram.py ligger direkt i subsystem_python.
    """
    return Path(__file__).resolve().parent


def default_results_root() -> Path:
    return project_root() / DEFAULT_RESULTS_DIRNAME


def default_system_dir(system: str) -> Path:
    return default_results_root() / system


def default_phase_curves_dir(system: str) -> Path:
    return default_system_dir(system) / DEFAULT_PHASE_CURVES_DIRNAME


def default_curves_path(system: str) -> Path:
    return default_phase_curves_dir(system) / DEFAULT_CURVES_NAME


def default_phase_diagram_dir(system: str) -> Path:
    return default_phase_curves_dir(system) / DEFAULT_PHASE_DIAGRAM_DIRNAME


def default_phase_diagram_path(system: str) -> Path:
    return default_phase_diagram_dir(system) / DEFAULT_PHASE_DIAGRAM_NAME


if __name__ == "__main__":
    main()