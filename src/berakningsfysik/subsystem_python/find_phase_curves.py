"""Denna fil hittar kompositioner som gäller vid binodal- och
spinodal-kurvorna, och skriver dessa till en csv-fil, curves.csv."""

from scipy.spatial import ConvexHull, convex_hull_plot_2d
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
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

    for index, T in enumerate(temps):
        comps_list = list(find_comps_at_temp(T, df, index))
        
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


def find_comps_at_temp(T, df, index):
    """Returnerar kompositioner for binodal
    och spinodal-kurvor vid temp T(som lista med floats)
    df = dataframe med data läst från 
    dataframe.csv
    """
    x_interpolated = df["x"]
    deltaG = df["deltaG_T"+str(index)]
    d2deltaG = df["d2deltaG_T"+str(index)]
    points = np.array(list(zip(x_interpolated, deltaG)))
    hull = ConvexHull(points)

    xa_list = []
    xb_list = []
    spinodal_xa_list = []
    spinodal_xb_list = []

    plt.clf()
    spinodals = []
    for simplex in hull.simplices:
        
        plt.plot(points[simplex, 0], points[simplex, 1], 'k-')

        start = simplex[0]  # ger index till startpunkten i simplex
        end = simplex[1]  # ger index till slutpunkten i simplex
        
        if start < end:
            x_start = points[start][0] #verkar vara samma som start?
            x_end = points[end][0]
            x_hull = np.array(x_interpolated[start: end + 1])
        
        else:
            x_start = points[end][0]
            x_end = points[start][0]
            x_hull = np.array(x_interpolated[end: start + 1])

        num_inflection_points = 0
        
        if x_start == 0 or x_start == 1:
            recent_sign = 1
        else:
            recent_sign = d2deltaG[min(start, end)-1]/abs(d2deltaG[min(start, end)-1])
       

        spinodal_comps = []

        for x_index, x in enumerate(x_hull):
            if d2deltaG[x_index+min(start, end)] < 0:
                if recent_sign == 1:
                    num_inflection_points += 1
                    spinodal_comps.append(x)
                    spinodals.append(x)
                recent_sign = -1

            else:
                if recent_sign == -1:
                    num_inflection_points += 1
                    spinodal_comps.append(x)
                    spinodals.append(x)
                recent_sign = 1
        
        if num_inflection_points == 2:

            deltaG_start = points[min(start,end)][1]
    
            if deltaG_start < 0:

                xa_list.append(float(x_start))
                xb_list.append(float(x_end))
                spinodal_xa_list.append(float(spinodals[0]))
                spinodal_xb_list.append(float(spinodals[1]))

    
    if not xa_list and np.any(deltaG > 0) and np.any(d2deltaG < 0):
    # gäller dessa villkor vet vi att ingen gemensam tangent finns
    # men spinodalkurvan finns ändå här 
    # detta gäller för låga temperaturer
        xa_list.append(0)
        xb_list.append(1)
        spinodal_xa_list.append(float(spinodals[0]))
        spinodal_xb_list.append(float(spinodals[1]))


    elif not xa_list:
    # gäller dessa villkor betyder det att ingen fasseparation sker
    # ingen binodal- eller spinodalkurva
    # i regel sker detta vid temperaturer över smältpunkten
        xa_list.append(np.nan)
        xb_list.append(np.nan)
        spinodal_xa_list.append(np.nan)
        spinodal_xb_list.append(np.nan)

    plt.savefig("plots/delta_G_mix/hull/hull_deltaG_T="+str(T))
    return xa_list, xb_list, spinodal_xa_list, spinodal_xb_list


if __name__ == "__main__":
    main()