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
    
    if n < 0:
        raise ValueError("Antalet mol kan inte vara negativt")

    if R != 8.314462618:
        raise ValueError("Gas konstanten kan inte vara något annat än 8.314462618 J/(mol*K)")

    entropy_sum = sum(x * math.log(x) for x in mole_fractions if x > 0)
    delta_s_mix = -R * n * entropy_sum
    return delta_s_mix

def calculate_delta_g_mix(mole_fractions, T, n, R):
    if T < 0:
        raise ValueError("Temperaturen kan inte vara negativ")
    
    delta_g_mix = - calculate_delta_s_mix(mole_fractions, n, R) * T
    return delta_g_mix

def calculate_delta_g_mix_with_h(mole_fractions, T, n, R, h):
    if T < 0:
        raise ValueError("Temperaturen kan inte vara negativ")
    
    delta_s_mix = calculate_delta_s_mix(mole_fractions, n, R)
    delta_g_mix_with_h = h - T * delta_s_mix
    return delta_g_mix_with_h