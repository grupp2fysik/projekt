"""Hjälpfunktioner för att hitta spinodalkurvan"""

#from thermodynamics import check_if_valid_x
from berakningsfysik.subsystem1.thermodynamics import check_if_valid_x

def entropy_second_derivative(x):
    check_if_valid_x(x)
    KB = 8.6173*10**(-5) # Boltzmanns konstant (eV/K)
    return -KB*1/(x*(1-x))