"""Denna fil använder datan i curves.csv för 
att plotta fasdiagrammet"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from build_dataframe import find_parameters
from find_phase_curves import columns

_, temps, alloy_name, qe_dir = find_parameters()

def main():

    print( "Plottar fasdiagram.")

    df = pd.read_csv("curves.csv")

    fig = plt.figure()

    plot_curve(df, "binodal", columns[1], columns[2])
    plot_curve(df, "spinodal", columns[3], columns[4])

    plt.xlim(0, 1)
    plt.ylim(0, temps[-1])
    plt.xlabel("x")
    plt.ylabel("T [K]")
    plt.title(f"Fasdiagram {alloy_name}")
    fig.legend()

    plt.savefig(f"plots/phase_diagrams/{alloy_name}_phase_diagram")


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



def turn_string_to_list(string):
    """Gör en sträng från curves.csv av formen
    "[2.897, 9.456]" till en lista med floats
    och returnerar denna. "nan" blir np.nan i 
    denna lista."""
    new_list = []
    comps_string = string.strip("[]")
    comps_list = comps_string.split(", ")
    for index, element in enumerate(comps_list):
        if element == "nan":
            comps_list[index] = np.nan
        else:
            comps_list[index] = float(element)
    
    return comps_list


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