"""Denna fil används endast för att generera spindodal- och
binodalgränser, för att testa plot_phase_diagram.py.
Kompositioenrna skrivs till generated_curves.csv"""

from scipy.spatial import ConvexHull, convex_hull_plot_2d
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from build_dataframe import temps

columns = ["T", "xa", "xb", "spinodal_xa", "spinodal_xb"]

def main():

    with open("generated_curves.csv", "w", newline="") as file:

        file.write(columns[0])

        for column in columns[1:]:
            file.write(","+column)
        file.write("\n")

        for T in temps:
            file.write(str(T)+",\n")

    file.close()
    xa_list = [[0, 0.37], [0, 0.4], [0, 0.42], [0.2, 0.5], [0.25, 0.55], [0.3, 0.6], [np.nan, np.nan]]
    xb_list = [[0.37, 1], [0.4, 1], [0.42, 1], [0.5, 0.9], [0.55, 0.7], [0.6, 0.8], [np.nan, np.nan]]
    spinodal_xa_list = [[0.1, 0.4], [0.15, 0.42], [0.2, 0.44], [0.3, 0.52], [0.35, 0.57], [0.4, 0.62], [np.nan, np.nan]]
    spinodal_xb_list = [[0.15, 0.45], [0.2, 0.47], [0.25, 0.49], [0.35, 0.57], [0.4, 0.62], [0.45, 0.67], [np.nan, np.nan]]

    df_curves = pd.read_csv("generated_curves.csv")

    df_curves["xa"] = xa_list
    df_curves["xb"] = xb_list
    df_curves["spinodal_xa"] = spinodal_xa_list
    df_curves["spinodal_xb"] = spinodal_xb_list

    df_curves.to_csv("generated_curves.csv", index=False)

if __name__ == "__main__":
    main()
