"""Denna fil läser från "dataframe.csv"
och plottar delta_S_mix för alla
temperaturer. Resultaten läggs i plots/delta_H_mix"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from build_dataframe import find_parameters, find_model

def print_deltaH_mix():

    print("Plottar entalpin.")

    _, _, alloy_name, qe_dir = find_parameters()
    _, x_values, y_values = find_model(qe_dir)
    df = pd.read_csv("dataframe.csv")
    x_interpolated = df["x"]
    figure = plt.figure()

    plt.clf()
    deltaH = df["deltaH"]
    points = np.array(list(zip(x_interpolated, deltaH)))
    plt.plot(points[:,0], points[:,1], 'k-')
    plt.plot(x_values, y_values, 'o', label="Mätpunkter från Quantum Espresso")

    plt.title("Redlich-Kister-polynom av \u0394H_mix")
    plt.xlabel("x")
    plt.ylabel("\u0394H_mix (eV/atom)")
    figure.legend()
    plt.savefig(f"plots/delta_H_mix/{alloy_name}_delta_H_mix")

if __name__ == "__main__":
    print_deltaH_mix()