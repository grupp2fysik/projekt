import math
R = 8.314462618
moles = 1.0

def calculate_delta_s_mix(mole_fractions, moles, R):
    if not math.isclose(sum(mole_fractions), 1.0, rel_tol=1e-99):
        raise ValueError("Summan av molbråken måste vara 1.0")
    
    for x in mole_fractions:
        if x < 0:
            raise ValueError("Molbråken kan inte vara negativ")

    entropy_sum = sum(x * math.log(x) for x in mole_fractions if x > 0)
    delta_s_mix = -R * moles * entropy_sum
    return delta_s_mix

def calculate_delta_g_mix(mole_fractions, T, moles, R):
    delta_g_mix = calculate_delta_s_mix(mole_fractions, moles, R) * T
    return delta_g_mix