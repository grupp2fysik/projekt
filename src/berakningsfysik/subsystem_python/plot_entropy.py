"""Denna fil läser från "dataframe.csv"
och plottar delta_S_mix. Resultaten läggs i plots/<legering>/delta_S_mix"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from build_dataframe import find_parameters
from parameters import *

def print_deltaS_mix():

    df = pd.read_csv(f"results/{alloy_name}/thermodynamics/dataframe.csv")
    x_interpolated = df["x"]
    figure = plt.figure()

    plt.clf()
    deltaS = df["deltaS"]
    points = np.array(list(zip(x_interpolated, deltaS)))
    plt.plot(points[:,0], points[:,1], 'k-')

    plt.title(f"\u0394S_mix ({alloy_name})")
    plt.xlabel("x")
    plt.ylabel("\u0394S_mix (eV/K per atom)")
    plt.savefig(f"plots/{alloy_name}/{alloy_name}_delta_S_mix")


if __name__ == "__main__":
    print_deltaS_mix()