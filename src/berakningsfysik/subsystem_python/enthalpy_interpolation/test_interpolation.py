import unittest
import math
import csv
import os

def check_negative_total_energy():
    """
    Kolla på om qe_outputs har negativa totala energin.
    """

    file_path0 = 'qe_outputs_test/qe_outputs0/enthalpy_dataset.csv'
    file_path1 = 'qe_outputs_test/qe_outputs1/enthalpy_dataset.csv'
    file_path2 = 'qe_outputs_test/qe_outputs2/enthalpy_dataset.csv'
    file_path3 = 'qe_outputs_test/qe_outputs3/enthalpy_dataset.csv'
    all_file_paths = [file_path0, file_path1, file_path2, file_path3]
    
    for index, number in enumerate(all_file_paths):
        with open(number, newline='') as csvfile:
            reader = csv.reader(csvfile, delimiter=',')
            header = next(reader)  # Hoppa över rubrikraden
            
            print("="*60)
            print(f"Kontrollerar negativa energier i '{header[2]}' för qe_outputs{index}.\n")
            
            negative_found = False
            for row in reader:
                x_value = row[1]  # x komposition
                energy = float(row[2])  # total_energy_Ry
                
                print(f"x={x_value}: {energy:.6f} Ry")
                
                if energy < 0:
                    negative_found = True
            
            if negative_found:
                print(f"\nNegativa energier detekterade (vilket förväntas för bundna system) i qe_outputs{index}.")
                print("="*60)
                print("\n")
            else:
                print("\nInga negativa energier finns.")
                print("="*60)

check_negative_total_energy()

def test_check_enthalpy_file():
    """
    Kolla på om filen existerar.
    """

    file_path0 = 'qe_outputs_test/qe_outputs0/enthalpy_dataset.csv'
    file_path1 = 'qe_outputs_test/qe_outputs1/enthalpy_dataset.csv'
    file_path2 = 'qe_outputs_test/qe_outputs2/enthalpy_dataset.csv'
    file_path3 = 'qe_outputs_test/qe_outputs3/enthalpy_dataset.csv'
    all_file_paths = [file_path0, file_path1, file_path2, file_path3]

    for index, number in enumerate(all_file_paths):
        if not os.path.isfile(number):
            raise FileNotFoundError(f"Filen enthalpy_dataset.csv hittades inte i qe_outputs{index}.")
        else:
            print(f"Filen enthalpy_dataset.csv finns och är tillgänglig i qe_outputs{index}.")

test_check_enthalpy_file()

def test_check_enthalpy_png():
    """
    Kontrollera att filen enthalpy_w_derivatives.png existerar.
    """

    file_path0 = 'qe_outputs_test/qe_outputs0/enthalpy_w_derivatives.png'
    file_path1 = 'qe_outputs_test/qe_outputs1/enthalpy_w_derivatives.png'
    file_path2 = 'qe_outputs_test/qe_outputs2/enthalpy_w_derivatives.png'
    file_path3 = 'qe_outputs_test/qe_outputs3/enthalpy_w_derivatives.png'
    all_file_paths = [file_path0, file_path1, file_path2, file_path3]

    for index, number in enumerate(all_file_paths):
        if not os.path.isfile(number):
            raise FileNotFoundError(f"Grafen enthalpy_w_derivatives.png hittades inte i qe_outputs{index}.")
        else:
            print(f"Grafen enthalpy_w_derivatives.png finns och är tillgänglig i qe_outputs{index}.")

test_check_enthalpy_png()

def test_check_directory():
    """
    Kolla på om katalogen finns.
    """

    directory_path0 = 'qe_outputs_test/qe_outputs0'
    directory_path1 = 'qe_outputs_test/qe_outputs1'
    directory_path2 = 'qe_outputs_test/qe_outputs2'
    directory_path3 = 'qe_outputs_test/qe_outputs3'
    all_directory_paths = [directory_path0, directory_path1, directory_path2, directory_path3]

    for number in all_directory_paths:
        if not os.path.isdir(number):
            raise NotADirectoryError(f"Katalogen '{number}' hittades inte.")
        else:
            print(f"Katalogen '{number}' finns och är tillgänglig.")

test_check_directory()

#def test_check_distance_of_square():



def test_check_L_coefficients():
    """
    Kolla på om L koefficienter existerar.
    Den kollar på om filen existerar för verifiering.
    """

    file_path = 'qe_outputs_test/qe_outputs3/rk_coeffs.npy'
    if not os.path.isfile(file_path):
        raise ValueError("L koefficienter finns inte")
    else:
        print("L koefficienter finns!")

test_check_L_coefficients()

# Testa att klarar av att vi skickar in konstigt fil (felmedellanden)
# Skicka ett directory som inte finns (felmedellanden)
# Tips på vad man kan göra istället.

# rot_min square error (hur långt det är)
# Se till att returnerar L koefficienter"""