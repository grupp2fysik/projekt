import numpy as np
from scipy.optimize import fsolve, minimize_scalar
from scipy.interpolate import interp1d
from dataclasses import dataclass
from typing import List, Tuple, Optional, Callable

@dataclass
class PhaseEquilibrium:
    R: 8.314462618
    T: float
    n: float
    x1: float
    x2: float
    G1: float
    G2: float

class GibbsMixing:
    def __init__(self, n: float, T: float, R = 8.314462618, omega: float):
        self.n = n
        self.T = T
        self.R = R
        self.omega = omega

    def calculate_ideal_mixing(self, x: float) -> float:
        if x < 0 or x > 1:
            raise ValueError("Composition x must be between 0 and 1")
        
        self.ideal = self.R * self.T * (x * np.log(x) + (1 - x) * np.log(1 - x)) / self.n
        return self.ideal
    
    def calculate_excess_mixing(self, x: float) -> float:
        if x < 0 or x > 1:
            raise ValueError("Composition x must be between 0 and 1")
        
        self.excess = self.omega * x * (1 - x)
        return self.excess
    
    def calculate_total_mixing(self, x: float) -> float:
        return self.ideal + self.excess
    
    def differentiate_total_mixing(self, x: float) -> float:
        d_ideal_dx = self.R * self.T * (np.log(x) - np.log(1 - x)) / self.n
        d_excess_dx = self.omega * (-2 * x + 1)
        return d_ideal_dx + d_excess_dx
    
    def differentiate_second_total_mixing(self, x: float) -> float:
        d2_ideal_dx2 = self.R * self.T * (1 / x + 1 / (1 - x)) / self.n
        d2_excess_dx2 = -2 * self.omega
        return d2_ideal_dx2 + d2_excess_dx2
    
class SpinodalCalculator:
    def __init__(self, gibbs_mixing: GibbsMixing):
        self.gibbs_mixing = gibbs_mixing

    def confirm_spinodal(self, x: float) -> bool:
        return self.gibbs_mixing.differentiate_second_total_mixing(x) == 0
    
    def find_spinodal(self, x_guess: float) -> Optional[float]:
        result = fsolve(self.gibbs_mixing.differentiate_second_total_mixing, x_guess)
        if self.confirm_spinodal(result[0]):
            return result[0]
        return None
    
class BinodalCalculator:
    def __init__(self, gibbs_mixing: GibbsMixing):
        self.gibbs_mixing = gibbs_mixing

        