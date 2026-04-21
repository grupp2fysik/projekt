"""Denna fil använder datan i curves.csv för 
att plotta fasdiagrammet"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from build_dataframe import temps

def main():
    df = pd.read_csv("curves.csv")

    fig = plt.figure()
    plt.xlim(0, 1)

    print(df.columns)

    print(df[df.columns[0]][0]) #printar 300 efterso detta är den första tempen
    print(df.index)

    for index in df.index:
        print(df[df.columns[0]][index])
        
    for column in df.columns:
        print(column)
        print("hej")
        
        for index in df[column].index:
            pass

    plt.savefig("plots/phase_diagram")
if __name__ == "__main__":
    main()