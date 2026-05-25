import argparse
import csv
import os
import pandas as pd
import numpy as np
from scipy.interpolate import interp1d
import subprocess
import sys

class PhaseDiagramAnalyser:
    def __init__(self, curves_file="results/TiAlN/phase_curves/curves.csv", config_file="alloy_parameters/TiAlN.csv"):
        """
        Initialisera PhaseDiagramAnalyser med data från en CSV-fil
        som innehåller information om binodala och spinodala kurvor.
        """

        self.curves_file = curves_file
        self.config_file = config_file
        self.data = pd.read_csv("results/TiAlN/phase_curves/curves.csv")
        self.df_curves = pd.read_csv(self.curves_file)
        self.df_config = pd.read_csv(self.config_file)

        self.xa_interpolated = interp1d(self.data["T"], self.data["xa"], kind="linear", fill_value="extrapolate")
        self.xb_interpolated = interp1d(self.data["T"], self.data["xb"], kind="linear", fill_value="extrapolate")

        self.spinodal_xa_interpolated = interp1d(self.data["T"], self.data["spinodal_xa"], kind="linear", fill_value="extrapolate")
        self.spinodal_xb_interpolated = interp1d(self.data["T"], self.data["spinodal_xb"], kind="linear", fill_value="extrapolate")

        self.specific_temperatures = []
        self.get_specific_temperatures()

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
            fraction_alpha = (x_beta - x) / (x_beta - x_alpha)
            fraction_beta = (x - x_alpha) / (x_beta - x_alpha)
            return fraction_alpha, fraction_beta

    def get_phase_compositions(self, T, x):
        """
        Bestäm faskompositionerna av α- och β-faserna
        baserat på den övergripande kompositionen x och temperaturen T.
        """
        x_alpha_binodal, x_beta_binodal = self.get_binodal_compositions(T)
        
        if np.isnan(x_alpha_binodal) or np.isnan(x_beta_binodal):
            return None, None

        if x <= x_alpha_binodal:
            return x, None
        elif x >= x_beta_binodal:
            return None, x
        else:
            return x_alpha_binodal, x_beta_binodal

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

        return [float(t) for t in self.data["T"].values]  # Convert to regular float
    
    def get_specific_temperatures(self, config_file="alloy_parameters/TiAlN.csv"):
        """
        Hämta specifika temperaturer från config-filen, om de finns.
        """
        specific_temperatures = []
        
        try:
            with open(config_file, 'r') as file:
                for line in file:
                    if line.startswith('specifika temperaturer,'):
                        temps_str = line.strip().split(',')[1]
                        temps_list = [float(temp.strip()) for temp in temps_str.split()]
                        specific_temperatures.extend(temps_list)
                        break      
            self.specific_temperatures = specific_temperatures
            
        except Exception as e:
            print(f"Warning: Could not read specific temperatures from {config_file}: {e}")
            self.specific_temperatures = []
            
    def analyze_temperature_range(self, start_temp, end_temp, compositions=None):
        """
        Analysera ett temperaturintervall från start_temp till end_temp.
        """

        if compositions is None:
            compositions = [0.1, 0.3, 0.5, 0.7, 0.9]

        all_temps = self.get_all_temperatures()
        available_temps = set(all_temps)
        single_temp_mode = (end_temp is None)
    
        if single_temp_mode:
            if self.specific_temperatures and start_temp not in self.specific_temperatures:
                print(f"{start_temp}K not in specific temperatures")
                print(f"Specific temperatures: {self.specific_temperatures}")
                return
        
            if start_temp not in available_temps:
                print(f"No temperature data available for {start_temp} K")
                print(f"Available temperatures: {sorted(available_temps)}")
                return
            
            end_temp = start_temp
        else:
            if start_temp > end_temp:
                print(f"Start temperature ({start_temp} K) is greater than end temperature ({end_temp} K). Please provide a valid temperature range.")
                return
    
        temps_in_range = [T for T in all_temps if start_temp <= T <= end_temp]
        
        for T in temps_in_range:
            print("="*60)
            print(f"\nPhase analysis at T = {T} K")
                
            for x in compositions:
                result = self.calculate_phase_properties(T, x)
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
        print("="*60)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Phase diagram analysis for materials')
    parser.add_argument('material', type=str, help='Material name (e.g., TiAlN)')
    parser.add_argument('temperatures', type=float, nargs='+', help='Temperature(s) in Kelvin (one or more)')
    parser.add_argument('--compositions', type=float, nargs='+', 
                       default=[0.1, 0.3, 0.5, 0.7, 0.9],
                       help='Compositions to analyze')
    
    args = parser.parse_args()
    curves_file = f"results/{args.material}/phase_curves/curves.csv"
    config_file = f"alloy_parameters/{args.material}.csv"
    
    try:
        analyzer = PhaseDiagramAnalyser(curves_file, config_file)
        
        for temp in args.temperatures:
            analyzer.analyze_temperature_range(temp, None, args.compositions)
            
    except Exception as e:
        print(f"Error: {e}")
        sys.exit(1)