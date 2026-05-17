import csv

import pandas as pd
import numpy as np
from scipy.interpolate import interp1d
import find_phase_curves_0
import subprocess
import sys

class PhaseDiagramAnalyser:
    def __init__(self, curves_file="results/TiAlN/phase_curves/curves.csv"):
        """
        Initialisera PhaseDiagramAnalyser med data från en CSV-fil
        som innehåller information om binodala och spinodala kurvor.
        """

        self.curves_file = curves_file
        self.data = pd.read_csv("results/TiAlN/phase_curves/curves.csv")
        self.df_curves = pd.read_csv(self.curves_file)

        self.xa_interpolated = interp1d(self.data["T"], self.data["xa"], kind="linear", fill_value="extrapolate")
        self.xb_interpolated = interp1d(self.data["T"], self.data["xb"], kind="linear", fill_value="extrapolate")

        self.spinodal_xa_interpolated = interp1d(self.data["T"], self.data["spinodal_xa"], kind="linear", fill_value="extrapolate")
        self.spinodal_xb_interpolated = interp1d(self.data["T"], self.data["spinodal_xb"], kind="linear", fill_value="extrapolate")

    def get_spinodal_compositions(self, T):
        """
        Bestäm kompositionerna av α- och β-faserna vid
        spinodala punkter baserat på temperaturen T.
        """

        if self.spinodal_xa_interpolated is None or self.spinodal_xb_interpolated is None:
            return None, None

        x_spinodal_alpha = self.spinodal_xa_interpolated(T)
        x_spinodal_beta = self.spinodal_xb_interpolated(T)
        return x_spinodal_alpha, x_spinodal_beta

    def get_binodal_compositions(self, T):
        """
        Bestäm kompositionerna av α- och β-faserna vid
        binodala punkter baserat på temperaturen T.
        """

        if self.xa_interpolated is None or self.xb_interpolated is None:
            return None, None

        x_alpha = self.xa_interpolated(T)
        x_beta = self.xb_interpolated(T)
        return x_alpha, x_beta

    def lever_rule(self, T, x):
        """
        Bestäm molfraktionerna eller massfraktionerna av α- och β-faserna baserat på den 
        övergripande kompositionen x och temperaturen T med hjälp av leverregeln.
        """

        x_alpha, x_beta = self.get_binodal_compositions(T)

        if np.isnan(x_alpha) or np.isnan(x_beta):
            return None, None

        if x <= x_alpha:
            return 1.0, 0.0
        elif x >= x_beta:
            return 0.0, 1.0
        else:        
            fraction_alpha = (x_beta - x) / (x_beta - x_alpha) #Cannot be zero since x is between x_alpha and x_beta or negative since x is less than x_beta
            fraction_beta = (x - x_alpha) / (x_beta - x_alpha) #Cannot be zero since x is between x_alpha and x_beta or negative since x is greater than x_alpha
            return fraction_alpha, fraction_beta
        
    def get_phase_compositions(self, T, x):
        """
        Bestäm faskompositionerna av α- och β-faserna
        baserat på den övergripande kompositionen x och temperaturen T.
        """

        x_alpha, x_beta = self.get_spinodal_compositions(T)

        if np.isnan(x_alpha) or np.isnan(x_beta):
            return None, None

        if x <= x_alpha:
            return x, None
        else:
            return x_alpha, x_beta

    def is_spinodal_decompositions(self, T, x):
        """
        Bestäm om systemet genomgår spinodal dekomposition baserat på temperatur T och komposition x.
        """

        x_spinodal_a, x_spinodal_b = self.get_spinodal_compositions(T)
        x_alpha, x_beta = self.get_binodal_compositions(T)

        if np.isnan(x_spinodal_a) or np.isnan(x_spinodal_b) or np.isnan(x_alpha) or np.isnan(x_beta):
            return False

        if x_spinodal_a <= x <= x_spinodal_b:
            return True
        elif x_alpha < x < x_spinodal_a or x_spinodal_b < x < x_beta:
            return False
        else:
            return False
    
    def is_binodal_decompositions(self, T, x):
        """
        Kolla om systemet genomgår binodal dekomposition baserat på temperatur T och komposition x.
        """

        x_alpha, x_beta = self.get_binodal_compositions(T)

        if np.isnan(x_alpha) or np.isnan(x_beta):
            return False

        if x < x_alpha or x > x_beta:
            return False
        else:
            return True
        
    def get_decomposition_type(self, T, x):
        """
        Bestäm vilken typ av dekomposition som sker 
        baserat på temperatur T och komposition x.
        """

        x_alpha, x_beta = self.get_binodal_compositions(T)
        x_spinodal_a, x_spinodal_b = self.get_spinodal_compositions(T)

        if np.isnan(x_alpha) or np.isnan(x_beta) or np.isnan(x_spinodal_a) or np.isnan(x_spinodal_b):
            return "Unknown (Insufficient data)"

        if x < x_alpha or x > x_beta:
            return "Single Phase (Stable)"
        elif x_spinodal_a <= x <= x_spinodal_b:
            return "Spinodal Decomposition (Unstable)"
        else:
            return "Binodal Decomposition (Metastable)"
        
    def calculate_phase_properties(self, T, x):
        """
        Bestäm fasegenskaperna baserat på temperatur T och komposition x,
        inklusive faskvoter, kompositioner och dekompositionsmekanism.
        """

        fraction_alpha, fraction_beta = self.lever_rule(T, x)
        composition_alpha, composition_beta = self.get_phase_compositions(T, x)
        mechanism = self.get_decomposition_type(T, x)
        return {
            'temperature': T,
            'overall_composition': x,
            'fraction_alpha': fraction_alpha,
            'fraction_beta': fraction_beta,
            'composition_alpha': composition_alpha,
            'composition_beta': composition_beta,
            'decomposition_mechanism': mechanism,
            'spinodal_decomposition': self.is_spinodal_decompositions(T, x)
        }

    def get_all_temperatures(self):
        """
        Hämta alla temperaturer från datan.
        """

        return self.data["T"].values
        
if __name__ == "__main__":
    data = pd.read_csv("results/TiAlN/phase_curves/curves.csv")
    analyzer = PhaseDiagramAnalyser("results/TiAlN/phase_curves/curves.csv")

    with open("alloy_parameters/TiAlN.csv", newline='') as csvfile:
        spamreader = csv.reader(csvfile, delimiter=',') #Läs in data från TiAlN.csv och lagra i en lista för att kunna uppdatera specifika temperaturer
        data = list(spamreader) #Lagra data i en lista för att kunna uppdatera specifika temperaturer

    atom_input = input("Skriva in atomnummer för element A (Ti): ")
    if not atom_input.strip().isdigit() or int(atom_input.strip()) <= 0:
        print(f"Ogiltig värde '{atom_input}'. Ange ett positivt heltal för atomnummer.")
        sys.exit(1) #Avsluta programmet om ogiltiga värden hittas

    min_input = input("Skriva in min temperatur i Kelvin K: ")
    if not min_input.strip().isdigit() or int(min_input.strip()) <= 0:
        print(f"Ogiltig värde '{min_input}'. Ange ett positivt heltal för min temperatur.")
        sys.exit(1)

    max_input = input("Skriva in max temperatur i Kelvin K: ")
    if not max_input.strip().isdigit() or int(max_input.strip()) <= 0:
        print(f"Ogiltig värde '{max_input}'. Ange ett positivt heltal för max temperatur.")
        sys.exit(1)

    temps_input = input("Skriva in temperaturer i Kelvin K (separerat med kommatecken): ")
    for temps in temps_input.split(","):
        if not temps.strip().isdigit():
            print(f"Ogiltig värde '{temps}'. Ange positiva heltal för temperaturer och/eller seperara med kommatecken.")
        if int(temps.strip()) < int(min_input.strip()) or int(temps.strip()) > int(max_input.strip()):
            print(f"Temperaturen måste vara mellan {min_input.strip()} K och {max_input.strip()} K. '{temps}' är utanför detta intervall.")
            sys.exit(1)
    
    temperature_list = []
    for temp in temps_input.split(","):
        if temp.strip() >= "0":
            temperature_list.append(temp.strip())
    new_temps = temps_input.split(",")

    for i, row in enumerate(data):
        if row and row[0].strip() == "atomer per metallplats":
            data[i][1] = atom_input.strip()
        if row and row[0].strip() == "specifika temperaturer":
            data[i][1] = " ".join(new_temps)
        if row and row[0].strip() == "temperatur max":
            data[i][1] = max_input.strip()
        if row and row[0].strip() == "temperatur min":
            data[i][1] = min_input.strip()
    
    with open("alloy_parameters/TiAlN.csv", 'w', newline='') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerows(data)

    print("="*60)
    print("Kör find_phase_curves_0.py för att generera fas kurvor...")
    print("="*60)

    try:
        result = subprocess.run([sys.executable, "find_phase_curves_0.py"], 
                              capture_output=True, 
                              text=True)
        
        if result.returncode == 0:
            print("✅ find_phase_curves_0.py kört framgångsrikt!")

        else:
            print(f"❌ Fel vid körning av find_phase_curves_0.py (Avslutningskod: {result.returncode})")
            print(f"Felutdata: {result.stderr}")
            sys.exit(1)
    except FileNotFoundError:
        print("❌ Error: find_phase_curves_0.py inte hittad.")
        print("Vänligen kontrollera att filen finns och är på rätt plats.")
        sys.exit(1)
    except Exception as e:
        print(f"❌ Oväntat fel: {e}")
        sys.exit(1)

    compositions = [0.1, 0.3, 0.5, 0.7, 0.9]
    analyzer = PhaseDiagramAnalyser("results/TiAlN/phase_curves/curves.csv")

    temperature_list_numeric = [int(float(temp)) for temp in temperature_list]

    for T in analyzer.get_all_temperatures():
        T_int = int(round(T))
        if T_int in temperature_list_numeric:
            print("="*60)
            print(f"\nPhase analysis at T = {T} K")
            
            for x in compositions:
                result = analyzer.calculate_phase_properties(T, x)
                print(f"\nOverall composition x = {x:.2f}")
                print(f"  Decomposition: {result['decomposition_mechanism']}")

                if result['fraction_alpha'] is not None:
                    print(f"  α phase fraction: {result['fraction_alpha']:.3f}")
                else:
                    print(f"  α phase fraction: None")    
                if result['fraction_beta'] is not None:
                    print(f"  β phase fraction: {result['fraction_beta']:.3f}")
                else:
                    print(f"  β phase fraction: None")
                
                if result['composition_alpha'] is not None:
                    print(f"  α phase composition: {result['composition_alpha']:.4f}")
                else:
                    print(f"  α phase composition: None")
                if result['composition_beta'] is not None:
                    print(f"  β phase composition: {result['composition_beta']:.4f}")
                else:
                    print(f"  β phase composition: None")
            print("\n")

    print("="*60)