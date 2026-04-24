import pandas as pd
import numpy as np
from scipy.interpolate import interp1d

class PhaseDiagramAnalyser:
    def __init__(self, curves_file="curves.csv"):
        self.curves_file = curves_file
        self.data = pd.read_csv("curves.csv")
        self.df_curves = pd.read_csv(self.curves_file)

        self.xa_interpolated = interp1d(self.data["T"], self.data["xa"], kind="linear", fill_value="extrapolate")
        self.xb_interpolated = interp1d(self.data["T"], self.data["xb"], kind="linear", fill_value="extrapolate")

        self.spinodal_xa_interpolated = interp1d(self.data["T"], self.data["spinodal_xa"], kind="linear", fill_value="extrapolate")
        self.spinodal_xb_interpolated = interp1d(self.data["T"], self.data["spinodal_xb"], kind="linear", fill_value="extrapolate")
        
    def get_binodal_comps(self, T):
        x_alpha = self.xa_interpolated(T)
        x_beta = self.xb_interpolated(T)
        return x_alpha, x_beta
    
    def get_spinodal_comps(self, T):
        x_spinodal_alpha = self.spinodal_xa_interpolated(T)
        x_spinodal_beta = self.spinodal_xb_interpolated(T)
        return x_spinodal_alpha, x_spinodal_beta
    
    def lever_rule(self, T, x):
        x_alpha, x_beta = self.get_binodal_comps(T)
        if x < x_alpha or x > x_beta:
            raise ValueError("Composition x is outside the binodal region.")
        
        fraction_alpha = (x_beta - x) / (x_beta - x_alpha) #Cannot be zero since x is between x_alpha and x_beta or negative since x is less than x_beta
        fraction_beta = (x - x_alpha) / (x_beta - x_alpha) #Cannot be zero since x is between x_alpha and x_beta or negative since x is greater than x_alpha
        
        return fraction_alpha, fraction_beta
    
    def find_phase_compositions(self, T, x):
        x_alpha, x_beta = self.get_binodal_comps(T)
        
        if x <= x_alpha:
            return x, None
        elif x >= x_beta:
            return None, x
        else:
            return x_alpha, x_beta
    
    def is_spinodal_compositions(self, T, x):
        x_spinodal_a, x_spinodal_b = self.get_spinodal_compositions(T)
        x_alpha, x_beta = self.get_binodal_compositions(T)

        if x_spinodal_a <= x <= x_spinodal_b:
            return True
        elif x_alpha < x < x_spinodal_a or x_spinodal_b < x < x_beta:
            return False
        else:
            return False
        
    