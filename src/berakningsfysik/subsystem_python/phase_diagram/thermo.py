"""
Denna fil innehåller funktioner för att beräkna termodynamiska egenskaper
relaterade till fasdiagram, såsom entropi av blandning och Gibbs fria energi av blandning.

Funktionerna tar hänsyn till molbråk, temperatur, antal mol och 
gas konstanten för att utföra beräkningarna. 

Det finns också en funktion som inkluderar entalpi av 
blandning i beräkningen av Gibbs fria energi.
"""

import math
R = 8.314462618
n = 1.0

def calculate_delta_s_mix(mole_fractions, n, R):
    """
    Beräkna entropi av blandning för en binär blandning.

    Parametrar:
    mole_fractions (list): En lista med molbråk för komponenterna i blandningen.
    n (float): Det totala antalet mol i blandningen.
    R (float): Gas konstanten (standard är 8.314462618 J/(mol*K)).
    """

    if not math.isclose(sum(mole_fractions), 1.0, rel_tol=1e-999):
        raise ValueError("Summan av molbråken måste vara 1.0")
    
    if type(mole_fractions) != list:
        raise TypeError("Molbråken måste vara en lista")
    
    for x in mole_fractions:
        if x < 0:
            raise ValueError("Molbråken kan inte vara negativ")
        if x > 1:
            raise ValueError("Molbråken kan inte vara större än 1")
    
    if n < 0:
        raise ValueError("Antalet mol kan inte vara negativt")

    if R != 8.314462618:
        raise ValueError("Gas konstanten kan inte vara något annat än 8.314462618 J/(mol*K)")

    entropy_sum = sum(x * math.log(x) for x in mole_fractions if 0 < x < 1)
    delta_s_mix = -R * n * entropy_sum
    return delta_s_mix

def calculate_delta_g_mix(mole_fractions, T, n, R):
    """
    Beräkna Gibbs fria energi av blandning för en binär blandning
    baserat på entropi av blandning och temperatur.

    Parametrar:
    mole_fractions (list): En lista med molbråk för komponenterna i blandningen.
    T (float): TTemperaturen i Kelvin.
    n (float): Det totala antalet mol i blandningen.
    R (float): Gas konstanten (standard är 8.314462618 J/(mol*K)).
    """

    if T < 0:
        raise ValueError("Temperaturen kan inte vara negativ")
    
    delta_g_mix = - calculate_delta_s_mix(mole_fractions, n, R) * T
    return delta_g_mix

def calculate_delta_g_mix_with_h(mole_fractions, T, n, R, delta_h):
    """
    Beräkna Gibbs fria energi av blandning för en binär blandning
    när entalpi av blandning är känd, baserat på entropi av blandning,
    temperatur och entalpi av blandning.

    Parametrar:
    mole_fractions (list): En lista med molbråk för komponenterna i blandningen.
    T (float): Temperaturen i Kelvin.
    n (float): Det totala antalet mol i blandningen.
    R (float): Gas konstanten (standard är 8.314462618 J/(mol*K)).
    delta_h (float): Entalpi av blandning i J/mol.
    """

    if T < 0:
        raise ValueError("Temperaturen kan inte vara negativ")
    
    delta_s_mix = calculate_delta_s_mix(mole_fractions, n, R)
    delta_g_mix_with_h = delta_h - T * delta_s_mix
    return delta_g_mix_with_h