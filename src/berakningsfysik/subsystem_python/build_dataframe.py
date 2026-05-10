"""Denna fil sparar interpolerad deltaH, deltaS och deltaG samt
andraderivatan av deltaG (för alla temp)
i en csv-fil som enkelt kan läsas till en pandas dataframe med kommando
pd.read_csv("dataframe.csv")"""

import math
import sys
import pandas as pd
import numpy as np
import argparse
from enthalpy_interpolation.interpolation import RedlichKisterModel, build_enthalpy_dataframe
from thermodynamics import entropy_per_atom
from spinodal_functions import entropy_second_derivative

alloy_name = sys.argv[1]
columns = ["x", "deltaH", "deltaS"]
num_of_inter_points = 500

def write_file():

    print(f"Räknar ut termodynamiska storheter för {alloy_name}. \nSkriver data till dataframe.csv")
 
    n, temps, _, qe_dir = find_parameters()
   
    x_interpolation = np.linspace(0.0, 1.0, num_of_inter_points)
    model, x_values, y_values = find_model(qe_dir)
    H_interpolated = find_deltaH(model, x_interpolation)

    # skriver alla kolumner samt interpolationspunkter
    with open("dataframe.csv", "w", newline="") as file:

        file.write(columns[0])

        for column in columns[1:]:
            file.write(","+column)

        for i in range(len(temps)):
            file.write(",deltaG_T"+str(i))

        file.write("\n")

        for x in x_interpolation:
            file.write(str(x)+",\n")

    file.close()

    write_data("deltaH", H_interpolated)
    write_data("deltaS", find_deltaS(2, x_interpolation))
    

    for i in range(len(temps)):
        write_data("deltaG_T"+str(i), find_deltaG(temps[i]))
        write_data("d2deltaG_T"+str(i), find_d2G(model, x_interpolation, temps[i], n))
        

def write_data(column, data_array):
    """Skriver datan i "data_array"
    under "column" i csv-filen.
    column = sträng från columns (t.ex "deltaH")
    data = vektor med data till kolumnen"""
    df = pd.read_csv("dataframe.csv")
    df[column] = data_array
    df.to_csv("dataframe.csv", index=False)

def find_deltaS(n, x_interpolated):
    """Hittar interpolerad deltaS_mix 
    n = atomer per metallplats"""
    deltaS = []
    for x in x_interpolated:
        deltaS.append(entropy_per_atom(x, n))
    return np.array(deltaS)

def find_deltaG(T):
    df = pd.read_csv("dataframe.csv")
    delta_S = df["deltaS"]
    delta_H = df["deltaH"]
    deltaG = delta_H - T*delta_S
    return deltaG

def find_model(qe_dir):
    parser = argparse.ArgumentParser(description="Skapa dataframe.")
    parser.add_argument(
        "data_dir",
        nargs="?",
        default="enthalpy_interpolation/qe_outputs/qe_outputs0",
        help="Katalog med .out-filer. Standard: qe_outputs",
    )
    parser.add_argument(
        "--glob",
        default="*x=*.out",
        help="Glob-mönster för filnamn. Standard: *x=*.out",
    )
    parser.add_argument(
        "--order",
        type=int,
        default=3,
        help="Ordning på Redlich-Kister-polynomet. Standard: 3",
    )

    args = parser.parse_args()

    #df = build_enthalpy_dataframe(args.data_dir, glob_pattern=args.glob)
    df = build_enthalpy_dataframe("enthalpy_interpolation/qe_outputs/" + qe_dir, "*x=*.out")

    x = df["x"].to_numpy()
    y = df["H_mix_eV_per_atom"].to_numpy()

    return RedlichKisterModel.fit(x=x, hmix=y, order=args.order), x, y

def find_deltaH(model, x_interpolation):
    """Hämtar interpolerade deltaH_mix
    x_interpolation = array med interpolationspunkter (använd np.linspace(...))"""
    
    return model.hmix(x_interpolation)

def find_d2G(model, x_interpolation, T, n):
    """returnerar vektor med 
    andraderivatan av deltaG"""
    d2deltaH = model.d2(x_interpolation)[1:-1]
    d2deltaS = []
    for x in x_interpolation[1:-1]:
        d2deltaS.append(entropy_second_derivative(x, n))

    d2deltaG = np.insert((d2deltaH - T*np.array(d2deltaS)), 0, np.nan)
    return np.append(d2deltaG, [np.nan])

def find_parameters():
    """Läser en csv-fil med materialspecifika 
    parametrar och returnerar dessa"""

    df = pd.read_csv(f"alloy_parameters/{alloy_name}.csv", header=None)
    parameters = df[1]
    n = int(parameters[0])
    T_start = int(parameters[1])
    T_end = int(parameters[2])
    qe_dir = parameters[3].strip()

    spec_temps = parameters[4]
    spec_temps_list = spec_temps.split()

    for index, temp in enumerate(spec_temps_list):
        spec_temps_list[index] = int(temp)
   
    
    temps = [i for i in range(T_start, T_end, math.ceil((T_end - T_start)/30))]
    
    temps += spec_temps_list
    temps.sort()
    return n, temps, alloy_name, qe_dir
    

if __name__ == "__main__":
    write_file()
