"""Denna fil beräknar delta_S_mix, delta_G_mix och delta_G_mix_with_h
för en blandning av två komponenter. Funktionerna tar hänsyn 
till molbråk, temperatur, antal mol och gas konstanten. 

Resultaten kan användas för att analysera termodynamiska 
egenskaper hos blandningen och för att förstå hur olika faktorer p
åverkar entropi och Gibbs fria energi i systemet."""

import math
R = 8.314462618
n = 1.0

def calculate_delta_s_mix(mole_fractions, n, R):
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
    if T < 0:
        raise ValueError("Temperaturen kan inte vara negativ")
    
    delta_g_mix = - calculate_delta_s_mix(mole_fractions, n, R) * T
    return delta_g_mix

def calculate_delta_g_mix_with_h(mole_fractions, T, n, R, delta_h):
    if T < 0:
        raise ValueError("Temperaturen kan inte vara negativ")
    
    delta_s_mix = calculate_delta_s_mix(mole_fractions, n, R)
    delta_g_mix_with_h = delta_h - T * delta_s_mix
    return delta_g_mix_with_h