"""Hjälpfunktioner för att hitta spinodalkurvan"""

def entropy_second_derivative(x):
    KB = 8.6173*10**(-5) # Boltzmanns konstant (eV/K)
    return -KB*1/(x*(1-x))