import pandas as pd
import numpy as np
from scipy.interpolate import interp1d
from build_dataframe import temps

def turn_string_to_list(string):
    """
    Convertera strängar av formen "[2.897, 9.456]" till en lista med floats.
    Hanterar även "nan" och "[nan]" genom att returnera np.nan i listan.

    Parameter:
    string (str): En sträng som representerar en lista av floats eller "nan".
    """

    if pd.isna(string) or string == "nan" or string == "[nan]":
        return [np.nan]
    
    comps_string = string.strip("[]")  #Ta bort hakparenteser
    if not comps_string:  #Om strängen är tom efter att ha tagit bort hakparenteser,
        return []
    
    comps_list = comps_string.split(", ")  #Dela strängen i en lista av element
    for index, element in enumerate(comps_list):
        if element == "nan":
            comps_list[index] = np.nan
        else:
            comps_list[index] = float(element)
    
    return comps_list

class PhaseDiagramAnalyser:
    def __init__(self, curves_file="curves.csv"):
        """
        Initialisera PhaseDiagramAnalyser med data från en CSV-fil
        som innehåller information om binodala och spinodala kurvor.

        Parameter:
        curves_file (str): Sökväg till CSV-filen som innehåller kurvdata.
        """

        self.curves_file = curves_file
        self.data = pd.read_csv("curves.csv")

        xa_lists = []
        xb_lists = []
        spinodal_xa_lists = []
        spinodal_xb_lists = []
        temperatures = []

        for index, row in self.data.iterrows():
            print(index)
            print(row)
            temperatures.append(row["T"])
            xa_lists.append(turn_string_to_list(row["xa"]))
            xb_lists.append(turn_string_to_list(row["xb"]))
            spinodal_xa_lists.append(turn_string_to_list(row["spinodal_xa"]))
            spinodal_xb_lists.append(turn_string_to_list(row["spinodal_xb"]))
        
        xa_values = [lst[0] if lst and not pd.isna(lst[0]) else np.nan for lst in xa_lists]
        xb_values = [lst[0] if lst and not pd.isna(lst[0]) else np.nan for lst in xb_lists]
        spinodal_xa_values = [lst[0] if lst and not pd.isna(lst[0]) else np.nan for lst in spinodal_xa_lists]
        spinodal_xb_values = [lst[0] if lst and not pd.isna(lst[0]) else np.nan for lst in spinodal_xb_lists]

        valid_indices = [i for i, t in enumerate(temperatures) 
                        if not np.isnan(xa_values[i]) and not np.isnan(xb_values[i]) and
                        not np.isnan(spinodal_xa_values[i]) and not np.isnan(spinodal_xb_values[i])]
        
        if valid_indices:
            valid_temps = [temperatures[i] for i in valid_indices]
            valid_xa = [xa_values[i] for i in valid_indices]
            valid_xb = [xb_values[i] for i in valid_indices]
            valid_spinodal_xa = [spinodal_xa_values[i] for i in valid_indices]
            valid_spinodal_xb = [spinodal_xb_values[i] for i in valid_indices]

            self.xa_interpolated = interp1d(valid_temps, valid_xa, kind="linear", fill_value="extrapolate")
            self.xb_interpolated = interp1d(valid_temps, valid_xb, kind="linear", fill_value="extrapolate")

            self.spinodal_xa_interpolated = interp1d(valid_temps, valid_spinodal_xa, kind="linear", fill_value="extrapolate")
            self.spinodal_xb_interpolated = interp1d(valid_temps, valid_spinodal_xb, kind="linear", fill_value="extrapolate")
        
        else:
            self.xa_interpolated = None
            self.xb_interpolated = None
            self.spinodal_xa_interpolated = None
            self.spinodal_xb_interpolated = None
            
    def get_spinodal_compositions(self, T):
        """
        Determine the compositions of the α and β phases 
        at the spinodal points based on temperature T.
        """

        if self.spinodal_xa_interpolated is None or self.spinodal_xb_interpolated is None:
            return np.nan, np.nan

        x_spinodal_alpha = self.spinodal_xa_interpolated(T)
        x_spinodal_beta = self.spinodal_xb_interpolated(T)
        return x_spinodal_alpha, x_spinodal_beta
    
    def get_binodal_compositions(self, T):
        """
        Determine the compositions of the α and β phases
        at the binodal points based on temperature T.
        """

        if self.xa_interpolated is None or self.xb_interpolated is None:
            return np.nan, np.nan

        x_alpha = self.xa_interpolated(T)
        x_beta = self.xb_interpolated(T)
        return x_alpha, x_beta

    def lever_rule(self, T, x):
        """
        Determine the mole fractions or the mass fractions of the α and β phases 
        based on the overall composition x and temperature T using the lever rule.
        """

        x_alpha, x_beta = self.get_binodal_compositions(T)

        if np.isnan(x_alpha) or np.isnan(x_beta):
            return np.nan, np.nan

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
        Determine the phase compositions of the α and β phases
        based on the overall composition x and temperature T.
        """

        x_alpha, x_beta = self.get_binodal_compositions(T)

        if np.isnan(x_alpha) or np.isnan(x_beta):
            return None, None
        
        if x <= x_alpha:
            return x, None
        elif x >= x_beta:
            return None, x
        else:
            return x_alpha, x_beta
    
    def is_spinodal_decompositions(self, T, x):
        """
        Determine if the system is undergoing spinodal decomposition based on temperature T and composition x.
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
        Check if the system is undergoing binodal decomposition based on temperature T and composition x.
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
        Determine the type of decomposition
        based on temperature T and composition x.
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
        Calculate the phase properties based on temperature T and composition x,
        including phase fractions, compositions, and decomposition mechanism.
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
    
if __name__ == "__main__":
    data = pd.read_csv("curves.csv")
    analyzer = PhaseDiagramAnalyser("curves.csv")
    
    compositions = [0.1, 0.3, 0.5, 0.7, 0.9]
    
    for T in temps:
        print("="*60)
        print(f"Phase analysis at T = {T} K")
        
        for x in compositions:
            result = analyzer.calculate_phase_properties(T, x)
            print(f"\nOverall composition x = {x:.2f}")
            print(f"  Decomposition: {result['decomposition_mechanism']}")
            
            if not np.isnan(result['fraction_alpha']):
                print(f"  α phase fraction: {result['fraction_alpha']:.3f}")
            else:
                print(f"  α phase fraction: None")
                
            if not np.isnan(result['fraction_beta']):
                print(f"  β phase fraction: {result['fraction_beta']:.3f}")
            else:
                print(f"  β phase fraction: None")
            
            if result['composition_alpha'] is not None and not np.isnan(result['composition_alpha']):
                print(f"  α phase composition: {result['composition_alpha']:.4f}")
            else:
                print(f"  α phase composition: None")
                
            if result['composition_beta'] is not None and not np.isnan(result['composition_beta']):
                print(f"  β phase composition: {result['composition_beta']:.4f}")
            else:
                print(f"  β phase composition: None")

    print("="*60)