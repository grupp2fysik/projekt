"""Denna fil hittar kompositioner som gäller vid binodal- och
spinodal-kurvorna, och skriver dessa till en csv-fil, curves.csv"""

from scipy.spatial import ConvexHull, convex_hull_plot_2d
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from build_dataframe import temps
from build_dataframe import temps

columns = ["T", "xa", "xb", "spinodal_xa", "spinodal_xb"]

def main():

    print("Letar efter binodal- och spinodalkurvor.")

    df = pd.read_csv("dataframe.csv")
    
    with open("curves.csv", "w", newline="") as file:

        file.write(columns[0])

        for column in columns[1:]:
            file.write(","+column)
        file.write("\n")

        for T in temps:
            file.write(str(T)+",\n")

    file.close()

    xa_list = []
    xb_list = []
    spinodal_xa_list = []
    spinodal_xb_list = []

    for T in temps:
        comps_list = list(find_comps_at_temp(T, df))
        
        xa_list.append(comps_list[0])
        xb_list.append(comps_list[1])
        spinodal_xa_list.append(comps_list[2])
        spinodal_xb_list.append(comps_list[3])

    df_curves = pd.read_csv("curves.csv")

    df_curves["xa"] = xa_list
    df_curves["xb"] = xb_list
    df_curves["spinodal_xa"] = spinodal_xa_list
    df_curves["spinodal_xb"] = spinodal_xb_list

    df_curves.to_csv("curves.csv", index=False)


def find_comps_at_temp(T, df):
    """Returnerar kompositioner for binodal
    och spinodal-kurvor vid temp T(som lista med np.float)
    df = dataframe med data läst från 
    dataframe.csv
    """

    x_interpolated = df["x"]
    deltaG = df["deltaG_T"+str(4)]
    d2deltaG = df["d2deltaG_T"+str(4)]
    points = np.array(list(zip(x_interpolated, deltaG)))
    hull = ConvexHull(points)

    xa_list = []
    xb_list = []
    spinodal_xa_list = []
    spinodal_xb_list = []

    for simplex in hull.simplices:

        start = simplex[0]  # ger index till startpunkten i simplex
        end = simplex[1]  # ger index till slutpunkten i simplex

        x_start = points[start][0]

        x_end = points[end][0]

        if end > start:
            x_hull = np.array(x_interpolated[start: end + 1])
            
        else:
            continue

        num_inflection_points = 0
        recent_sign = 1

        spinodal_comps = []
        for index, x in enumerate(x_hull):

            if d2deltaG[index+start] < 0:
                if recent_sign == 1:
                    num_inflection_points += 1
                    spinodal_comps.append(x)
                recent_sign = -1

            else:
                if recent_sign == -1:
                    num_inflection_points += 1
                    spinodal_comps.append(x)
                recent_sign = 1
        
        if num_inflection_points == 2:
            xa_list.append(x_start)
            xb_list.append(x_end)
            spinodal_xa_list.append(spinodal_comps[0])
            spinodal_xb_list.append(spinodal_comps[1])
        
    return xa_list, xb_list, spinodal_xa_list, spinodal_xb_list



if __name__ == "__main__":
    main()