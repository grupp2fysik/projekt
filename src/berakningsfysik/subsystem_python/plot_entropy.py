"""Denna fil läser från "dataframe.csv"
och plottar delta_S_mix för alla
temperaturer. Resultaten läggs i delta_S_mix"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from build_dataframe import find_parameters
from parameters import *

def print_deltaS_mix():

    print("Plottar entropin.")

    #_, _, alloy_name, qe_dir = find_parameters()
    df = pd.read_csv("dataframe.csv")
    x_interpolated = df["x"]
    figure = plt.figure()

    plt.clf()
    deltaS = df["deltaS"]
    points = np.array(list(zip(x_interpolated, deltaS)))
    plt.plot(points[:,0], points[:,1], 'k-')

    plt.title(f"\u0394S_mix ({alloy_name})")
    plt.xlabel("x")
    plt.ylabel("\u0394S_mix (eV/K per atom)")
    plt.savefig(f"plots/delta_S_mix/{alloy_name}_delta_S_mix")

def return_points(filename):
    df = pd.read_csv(filename)
    x_interpolated = df["x"]

if __name__ == "__main__":
    print_deltaS_mix()