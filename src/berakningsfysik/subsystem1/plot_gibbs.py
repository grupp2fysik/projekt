"""Denna fil läser från "dataframe.csv"
och plottar delta_G_mix för alla
temperaturer. Resultaten läggs i plots/delta_G_mix"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.spatial import ConvexHull, convex_hull_plot_2d
from build_dataframe import temps


def print_deltaG_mix():

    df = pd.read_csv("dataframe.csv")

    x_interpolated = df["x"]
    figure = plt.figure()

    for index, T in enumerate(temps):

        plt.clf()
        deltaG = df["deltaG_T"+str(index)]
        d2deltaG = df["d2deltaG_T"+str(index)]

        points = np.array(list(zip(x_interpolated, deltaG)))
        d2points = np.array(list(zip(x_interpolated, d2deltaG)))

        hull = ConvexHull(points)

        plt.plot(points[:,0], points[:,1], 'k-')
       
    #for simplex in hull.simplices:
      #  plt.plot(points[simplex, 0], points[simplex, 1], 'k-') # hämta kolumn 0/1 för värdena som ingår i simplex

        plt.title("delta_G_mix för T = "+str(T)+"K")
        plt.xlabel("x")
        plt.ylabel("delta_G_mix")
        plt.savefig("plots/delta_G_mix/delta_G_mix_T="+str(T)+"K")

        plt.clf()

        plt.plot(d2points[:,0], d2points[:,1], 'b-')
        plt.plot(d2points[:,0], np.zeros(d2points[:,1].size), 'k--')
        plt.title("andraderivata delta_G_mix för T = "+str(T)+"K")
        plt.xlabel("x")
        plt.ylabel("andraderivata av delta_G_mix")
        plt.savefig("plots/delta_G_mix/d2delta_G_mix_T="+str(T)+"K")


if __name__ == "__main__":
    print_deltaG_mix()