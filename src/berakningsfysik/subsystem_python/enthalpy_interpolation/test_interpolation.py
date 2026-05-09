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

check_negative_total_energy()

#def test_check_file():
    


"""class InterpolationTests(unittest.TestCase):
    En test för...

    #def test_check_file(self):

    def check_energy_negative(self):
        with open('enthalpy_dataset.csv', mode='r') as file: 
            # Create a CSV reader object 
            csv_reader = csv.reader(file)

            # Skip the header row (if there is one) 
            next(csv_reader, None) 

            for row in csv_reader: 
                print(row) 


    #def test_check_directory(self):


    #def test_check_distance_of_square(self):


    #def test_check_L_coefficients(self):

        

# Testa att klarar av att vi skickar in konstigt fil (felmedellanden)
# Skicka ett directory som inte finns (felmedellanden)
# Tips på vad man kan göra istället.

# rot_min square error (hur långt det är)
# Se till att returnerar L koefficienter"""