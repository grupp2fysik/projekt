"""Denna fil använder datan i curves.csv för 
att plotta fasdiagrammet"""

import pandas as pd
import numpy as np
from pathlib import Path
import matplotlib.pyplot as plt
from build_dataframe import temps
from find_phase_curves_0 import columns

DEFAULT_SYSTEM = "TiAlN"
DEFAULT_RESULTS_DIRNAME = "results"
DEFAULT_PHASE_CURVES_DIRNAME = "phase_curves"
DEFAULT_PHASE_DIAGRAM_DIRNAME = "phase_diagram"
DEFAULT_CURVES_NAME = "curves.csv"
DEFAULT_PHASE_DIAGRAM_NAME = "phase_diagram.png"


def project_root() -> Path:
    """
    Returnerar subsystem_python-mappen.
    plot_phase_diagram_0.py ligger direkt i subsystem_python.
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

    
def main():

    print("Plottar fasdiagram.")

    curves_path = default_curves_path(DEFAULT_SYSTEM)
    output_path = default_phase_diagram_path(DEFAULT_SYSTEM)

    output_path.parent.mkdir(parents=True, exist_ok=True)

    print(f"Läser kurvor från: {curves_path}")
    print(f"Sparar fasdiagram till: {output_path}")

    df = pd.read_csv(curves_path)

    fig = plt.figure()

    plot_curve(df, "binodal", columns[1], columns[2])
    plot_curve(df, "spinodal", columns[3], columns[4])

    plt.xlim(0, 1)
    plt.ylim(0,)
    plt.xlabel("x")
    plt.ylabel("T [K]")
    plt.title("Fasdiagram")
    fig.legend()

    plt.savefig(output_path, dpi=200)


def plot_curve(df, curve_type, column1, column2):
    """Sätter ihop kolumnerna column1 och column2 
    (t.ex xa och xb)
    och plottar mot temperaturen.
    df är dataframe från curves.csv
    curve_type är "spinodal" eller "binodal" """

    line_styles = {"spinodal": "r--", "binodal": "b-"}

    main_list1, max_num_elements1 = return_comps_main_list(df, column1)
    main_list2, max_num_elements2 = return_comps_main_list(df, column2)

    main_list2.reverse()
    max_num_elements = max(max_num_elements1, max_num_elements2)

    for j in range(0, max_num_elements):

            x_values = []
            y_values =[]

            for index, comps_list in enumerate(main_list1):
                if len(comps_list) > j:
                    if not np.isnan(comps_list[j]):
                        x_values.append(comps_list[j])
                        y_values.append(temps[index])
                
            for index, comps_list in enumerate(main_list2):
                if len(comps_list) > j:
                    if not np.isnan(comps_list[j]):
                        x_values.append(comps_list[j])
                        y_values.append(temps[-index-1])
                
            plt.plot(x_values, y_values, line_styles[curve_type], label = curve_type)


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
    """Tar en kolumn från curves.csv, t.ex "xa",
    och omvandlar den från sträng till en nästlad lista
    (som ser ut som kolumnen i csv-filen). 
    Returnerar listan och längden på den av de inre
    listorna som är längst."""

    comps_main_list = []
    max_num_elements = 1
    for index in df.index:

        comps_list = turn_string_to_list(df[column][index])
        comps_main_list.append(comps_list)
        if len(comps_list) > max_num_elements:
            max_num_elements = len(comps_list)

    return comps_main_list, max_num_elements


if __name__ == "__main__":
    main()