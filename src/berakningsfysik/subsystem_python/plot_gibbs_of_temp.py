"""
Denna fil läser från "dataframe.csv"
och plottar delta_G_mix för alla
temperaturer som givits som argument i terminalen.
Resultaten läggs i plots/<legering>/delta_G_mix.
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from build_dataframe import find_parameters
import sys
from thermodynamics import check_if_T_in_temps, T_string, check_if_valid_T
from parameters import *

def print_delta_G_mix():
    """
    Plottar delta_G_mix för alla temperaturer som givits som argument i terminalen.
    """

    df = pd.read_csv(f"results/{alloy_name}/thermodynamics/dataframe.csv")

    spec_temps = sys.argv[2:]

    if not spec_temps:
        raise Exception("Ange vilken temperatur du vill analysera när du anropar filen.\
\nSe README.txt för instruktioner.")

    x_interpolated = df["x"]
    figure = plt.figure()

    for T in spec_temps:
        
        check_if_valid_T(T)
        T = float(T)
        print(f"Plottar Gibbs fria blandningsenergi för T = {T}K")

        check_if_T_in_temps(T, temps)
        plt.clf()
        
        index = temps.index(T)
        
        deltaG = df["deltaG_T"+str(index)]
        
        points = np.array(list(zip(x_interpolated, deltaG)))
        
        plt.plot(points[:,0], points[:,1], 'k-')
        
        plt.title("\u0394G_mix för T = "+str(T)+f"K ({alloy_name})")
        plt.xlabel("x")
        plt.ylabel("\u0394G_mix (eV/atom)")
        plt.savefig(f"plots/{alloy_name}/delta_G_mix/{alloy_name}_delta_G_mix_T="+T_string(T)+"K")

        plt.clf()

if __name__ == "__main__":
    print_delta_G_mix()