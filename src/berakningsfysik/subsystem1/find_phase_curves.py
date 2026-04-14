"""Denna fil hittar kompositioner som gäller vid binodal- och
spinodal-kurvorna, och skriver dessa till en csv-fil"""

from scipy.spatial import ConvexHull, convex_hull_plot_2d
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from build_dataframe import temps

def main():

    df = pd.read_csv("dataframe.csv")
    x_interpolated = df["x"]
    deltaG = df["deltaG_T"+str(4)]
    d2deltaG = df["d2deltaG_T"+str(4)]
    points = np.array(list(zip(x_interpolated, deltaG)))

    hull = ConvexHull(points)
    print("längd d2G = ", d2deltaG.size)

    for simplex in hull.simplices:

        plt.plot(points[simplex, 0], points[simplex, 1], 'k-')

        start = simplex[0]  # ger index till startpunkten
        end = simplex[1]  # ger index till slutpunkten

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
                    print("inflektionspunkt: ", x)
                    num_inflection_points += 1
                    spinodal_comps.append(x)
                recent_sign = 1
        
        if num_inflection_points == 2:
            print("xa_s = ", spinodal_comps[0], "\nxb_s = ", spinodal_comps[1])   
            print("xa = ", x_start, "xb =", x_end) 

    plt.savefig("convex_hull")       



if __name__ == "__main__":
    main()