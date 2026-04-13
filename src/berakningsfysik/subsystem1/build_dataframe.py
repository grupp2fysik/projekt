import pandas as pd
import numpy as np
import argparse
from enthalpy_interpolation.interpolation import RedlichKisterModel, build_enthalpy_dataframe
from thermodynamics import entropy_per_atom

columns = ["Index", "deltaH", "deltaS"]
temps = [50, 150, 250, 350]

def write_file():
    parser = argparse.ArgumentParser(description="Skapa dataframe.")
    parser.add_argument(
        "data_dir",
        nargs="?",
        default="enthalpy_interpolation/qe_outputs",
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

    df = build_enthalpy_dataframe(args.data_dir, glob_pattern=args.glob)

    x = df["x"].to_numpy()
    y = df["H_mix_eV_per_atom"].to_numpy()

    model = RedlichKisterModel.fit(x=x, hmix=y, order=args.order)

    x_grid = np.linspace(0.0, 1.0, 5)
    H_grid = model.hmix(x_grid)

    
    # skapa alla kolumner samt index/interpolationspunkter
    with open("dataframe.csv", "w", newline="") as file:
        file.write(columns[0])
        for column in columns[1:]:
            file.write(","+column)
        for i in range(len(temps)):
            file.write(",deltaG_T"+str(i))
        file.write("\n")
        for x in x_grid:
            file.write(str(x)+",\n")
    file.close()

    write_data("deltaH", H_grid)
    write_data("deltaS", find_deltaS(2, x_grid))

    for i in range(len(temps)):
        write_data("deltaG_T"+str(i), find_deltaG(temps[i]))
        

def write_data(column, data_array):
    """column = sträng från columns (t.ex deltaH)
    data = vektor med data till kolumnen"""
    df = pd.read_csv("dataframe.csv")
    print(df)
    df[column] = data_array
    print(df)
    df.to_csv("dataframe.csv", index=False)

def find_deltaS(n, x_grid):
    deltaS = []
    for x in x_grid:
        deltaS.append(entropy_per_atom(x, n))
    print(deltaS)
    return np.array(deltaS)

def find_deltaG(T):
    df = pd.read_csv("dataframe.csv")
    delta_S = df["deltaS"]
    delta_H = df["deltaH"]
    deltaG = delta_H - T*delta_S
    return deltaG

if __name__ == "__main__":
    write_file()
