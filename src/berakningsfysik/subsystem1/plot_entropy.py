"""Denna fil läser från "dataframe.csv"
och plottar delta_S_mix för alla
temperaturer. Resultaten läggs i delta_S_mix"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from build_dataframe import temps

def print_deltaS_mix():

    df = pd.read_csv("dataframe.csv")
    x_interpolated = df["x"]
    figure = plt.figure()

    plt.clf()
    deltaS = df["deltaS"]
    points = np.array(list(zip(x_interpolated, deltaS)))
    plt.plot(points[:,0], points[:,1], 'k-')
    plt.plot(points[:,0], np.zeros(points[:,1].size), 'k--')

    plt.title("delta_S_mix (oberoende av temperatur)")
    plt.xlabel("x")
    plt.ylabel("delta_S_mix")
    plt.savefig("delta_S_mix/delta_S_mix")

def return_points(filename):
    df = pd.read_csv(filename)
    x_interpolated = df["x"]

if __name__ == "__main__":
    print_deltaS_mix()