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

    print("kolumner: ", df.columns)

    print(df[df.columns[0]][0]) #printar 300 eftersom detta är den första tempen
    print("index: ", df.index)

    for index in df.index:
        
        temp = float(df[df.columns[0]][index])
        print("temp: ", temp)
        
        #plt.plot(np.arange([0, 0.2, 0.3, 0.5, 0.6, 0.9]), temps)
        for column in df.columns[1:]:
            comps = df[column][index]
            new_comps = turn_string_to_list(comps)
            
            #print(comps.size)
            plt.plot(new_comps, [temp]*len(new_comps), "o")
            
    for column in df.columns[1:]:
        comps_dataframe = df[column]
        comps = []
        max_num_comps = 1
        for index, element in enumerate(comps_dataframe):
            element_list = turn_string_to_list(element)
            for comps in element_list:
                
            if len(element_list) > max_num_comps:
                max_num_comps = len(element_list)
            comps.append(element_list)
            print(element_list)
        


            

        
        

    plt.savefig("plots/phase_diagram")

def turn_string_to_list(string):
    """Gör en sträng från curves.csv av formen
    "[2.897, 9.456]" till en lista med floats
    och returnerar denna. "nan" blir np.nan i 
    denna lista."""
    new_list = []
    comps_string = string.strip("[]")
    comps_list = comps_string.split(", ")
    for index, element in enumerate(comps_list):
        if element == "nan":
            comps_list[index] = np.nan
        else:
            comps_list[index] = float(element)
    
    return comps_list

if __name__ == "__main__":
    main()