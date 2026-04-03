"""Hjälpfunktioner för att hitta spinodalkurvan"""

#from thermodynamics import check_if_valid_x
from berakningsfysik.subsystem1.thermodynamics import check_if_valid_x
from berakningsfysik.subsystem1.thermodynamics import entropy_at_endpoint

def entropy_second_derivative(x):
    check_if_valid_x(x)
    if entropy_at_endpoint == 0:
        raise Exception("Derivatan inte definierad vid ändpunkterna.")
    KB = 8.6173*10**(-5) # Boltzmanns konstant (eV/K)
    return -KB*1/(x*(1-x))