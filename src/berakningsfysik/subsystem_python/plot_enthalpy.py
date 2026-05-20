"""
Denna fil läser från "dataframe.csv"
och plottar delta_H_mix för alla
temperaturer. Resultaten läggs i plots/legering/delta_H_mix
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from build_dataframe import find_parameters, find_model
from parameters import *

def print_delta_H_mix():
    """
    Plottar delta_H_mix för alla temperaturer. Resultaten läggs i plots/legering/delta_H_mix
    """
    
    print("Plottar entalpin.")

    qe_calc_df = pd.read_csv(qe_dir)
    df = pd.read_csv(f"results/{alloy_name}/thermodynamics/dataframe.csv")
    x_interpolated = df["x"]
    figure = plt.figure()

    plt.clf()
    deltaH = df["deltaH"]
    points = np.array(list(zip(x_interpolated, deltaH)))

    x_values = qe_calc_df["x"]
    y_values = qe_calc_df["hmix_ev_per_atom"]

    plt.plot(points[:,0], points[:,1], 'k-')
    plt.plot(x_values, y_values, 'o', label="Mätpunkter från Quantum Espresso")

    plt.title(f"Redlich-Kister-polynom av \u0394H_mix ({alloy_name})")
    plt.xlabel("x")
    plt.ylabel("\u0394H_mix (eV/atom)")
    figure.legend()
    plt.savefig(f"plots/{alloy_name}/{alloy_name}_delta_H_mix")

if __name__ == "__main__":
    print_delta_H_mix()