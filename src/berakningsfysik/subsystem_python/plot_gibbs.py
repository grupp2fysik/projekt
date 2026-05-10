"""Denna fil läser från "dataframe.csv"
och plottar delta_G_mix för alla
temperaturer. Resultaten läggs i plots/delta_G_mix/<material</delta_G_mix"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from build_dataframe import find_parameters

def print_deltaG_mix():

    parameters = find_parameters()
    df = pd.read_csv("dataframe.csv")

    x_interpolated = df["x"]
    figure = plt.figure()

    for index, T in enumerate(temps):

        plt.clf()
        deltaG = df["deltaG_T"+str(index)]
        d2deltaG = df["d2deltaG_T"+str(index)]

        points = np.array(list(zip(x_interpolated, deltaG)))
        d2points = np.array(list(zip(x_interpolated, d2deltaG)))

        plt.plot(points[:,0], points[:,1], 'k-')
        
        plt.title("delta_G_mix för T = "+str(T)+"K")
        plt.xlabel("x")
        plt.ylabel("delta_G_mix")
        plt.savefig(f"plots/delta_G_mix/{alloy_name}/delta_G_mix/delta_G_mix_T="+str(T)+"K")

        plt.clf()

        plt.plot(d2points[:,0], d2points[:,1], 'b-')
        plt.plot(d2points[:,0], np.zeros(d2points[:,1].size), 'k--')
        plt.title("andraderivata delta_G_mix för T = "+str(T)+"K")
        plt.xlabel("x")
        plt.ylabel("andraderivata av delta_G_mix")
        plt.savefig(f"plots/delta_G_mix/{alloy_name}/second_derivative/d2delta_G_mix_T="+str(T)+"K")

def return_points(filename):
    df = pd.read_csv(filename)
    x_interpolated = df["x"]

if __name__ == "__main__":
    print_deltaG_mix()