import unittest
import math
import csv
import os

def count_number_of_catalogues_in_directory(directory_path):
    """
    Räkna antalet kataloger i en katalog.
    """
    
    if not os.path.isdir(directory_path):
        raise NotADirectoryError(f"Katalogen '{directory_path}' hittades inte.")
    
    count = 0
    for entry in os.listdir(directory_path):
        entry_path = os.path.join(directory_path, entry)
        if os.path.isdir(entry_path):
            count += 1
    return count
    
def verify_directories_exist(directory_path):
    count = count_number_of_catalogues_in_directory(directory_path)

    all_directory_paths = []
    for i in range(count):
        file_path = f'qe_outputs_test/qe_outputs{i}'
        all_directory_paths.append(file_path)

    for number in all_directory_paths:
        if not os.path.isdir(number):
            raise NotADirectoryError(f"Katalogen '{number}' hittades inte.")
        else:
            print(f"Katalogen '{number}' finns och är tillgänglig.")
    
    return count

count = count_number_of_catalogues_in_directory('qe_outputs_test')
verify_directories_exist('qe_outputs_test')

def check_negative_total_energy():
    """
    Kolla på om qe_outputs har negativa totala energin.
    """

    all_file_paths = []
    for i in range(count):
        file_path = f'qe_outputs_test/qe_outputs{i}/enthalpy_dataset.csv'
        all_file_paths.append(file_path)
        
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

    all_file_paths = []
    for i in range(count):
        file_path = f'qe_outputs_test/qe_outputs{i}/enthalpy_dataset.csv'
        all_file_paths.append(file_path)

    for index, number in enumerate(all_file_paths):
        if not os.path.isfile(number):
            raise FileNotFoundError(f"Filen enthalpy_dataset.csv hittades inte i qe_outputs{index}.")
        else:
            print(f"Filen enthalpy_dataset.csv finns och är tillgänglig i qe_outputs{index}.")

def test_check_enthalpy_png():
    """
    Kontrollera att filen enthalpy_w_derivatives.png existerar
    i alla qe_outputs mappar.
    """

    all_file_paths = []
    for i in range(count):
        file_path = f'qe_outputs_test/qe_outputs{i}/enthalpy_w_derivatives.png'
        all_file_paths.append(file_path)

    for index, number in enumerate(all_file_paths):
        if not os.path.isfile(number):
            raise FileNotFoundError(f"Grafen enthalpy_w_derivatives.png hittades inte i qe_outputs{index}.")
        else:
            print(f"Grafen enthalpy_w_derivatives.png finns och är tillgänglig i qe_outputs{index}.")

    #def test_check_distance_of_square():
test_check_enthalpy_png()

def test_check_L_coefficients():
    """
    Kolla på om L koefficienter existerar.
    Den kollar på om filerna existerar för verifiering.
    """

    all_file_paths = []
    for i in range(count):
        directory_path = f'qe_outputs_test/qe_outputs{i}/rk_coeffs.npy'
        all_file_paths.append(directory_path)

    for index, number in enumerate(all_file_paths):
        if not os.path.isfile(number):
            raise ValueError(f"Filen rk_coeffs.npy finns inte i qe_outputs{index}, därför finns det inga L koefficienter.")
        else:
            print(f"Filen rk_coeffs.npy finns i qe_outputs{index}, därför finns det L koefficienter!")

test_check_L_coefficients()

# Testa att klarar av att vi skickar in konstigt fil (felmedellanden)
# Skicka ett directory som inte finns (felmedellanden)
# Tips på vad man kan göra istället.

# rot_min square error (hur långt det är)
# Se till att returnerar L koefficienter"""