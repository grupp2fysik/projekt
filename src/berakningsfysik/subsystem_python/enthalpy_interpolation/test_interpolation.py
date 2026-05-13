import unittest
import math
import csv
import os

def check_negative_total_energy():
    with open('qe_outputs/qe_outputs3/enthalpy_dataset.csv', newline='') as csvfile:
        reader = csv.reader(csvfile, delimiter=',')
        header = next(reader)  # Hoppa över rubrikraden
        
        print(f"Checking for negative energies in '{header[2]}'.\n")
        
        negative_found = False
        for row in reader:
            x_value = row[1]  # x komposition
            energy = float(row[2])  # total_energy_Ry
            
            print(f"x={x_value}: {energy:.6f} Ry")
            
            if energy < 0:
                negative_found = True
        
        if negative_found:
            print("\nNegativa energier detekterade (vilket förväntas för bundna system)")
        else:
            raise ValueError("\nInga negativa energier finns.")

def test_check_file():
    file_path = 'qe_outputs/qe_outputs3/enthalpy_dataset.csv'
    if not os.path.isfile(file_path):
        raise FileNotFoundError(f"Filen '{file_path}' hittades inte.")
    else:
        print(f"Filen '{file_path}' finns och är tillgänglig.")

def test_check_directory():
    directory_path = 'qe_outputs/qe_outputs3'
    if not os.path.isdir(directory_path):
        raise NotADirectoryError(f"Katalogen '{directory_path}' hittades inte.")
    else:
        print(f"Katalogen '{directory_path}' finns och är tillgänglig.")

def test_check_distance_of_square():
    


#def test_check_L_coefficients():


        

# Testa att klarar av att vi skickar in konstigt fil (felmedellanden)
# Skicka ett directory som inte finns (felmedellanden)
# Tips på vad man kan göra istället.

# rot_min square error (hur långt det är)
# Se till att returnerar L koefficienter"""