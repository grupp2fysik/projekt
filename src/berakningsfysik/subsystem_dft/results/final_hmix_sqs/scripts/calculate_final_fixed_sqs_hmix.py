RY_TO_EV = 13.605693122994

# Referensenergier från selected-beräkningarna
E_TiN_total = -848.67899733
E_AlN_total = -271.61076090

atoms_TiN = 8
atoms_AlN = 8
atoms_TiAlN = 64

E_TiN_atom = E_TiN_total / atoms_TiN
E_AlN_atom = E_AlN_total / atoms_AlN

data = [
    (0.125, -6212.17689898),
    (0.25,  -5634.99155838),
    (0.375, -5057.77759277),
    (0.5,   -4480.66132264),
    (0.625, -3903.57203596),
    (0.75,  -3326.55831569),
    (0.875, -2749.59994982),
]

print("x,total_energy_ry,hmix_ry_per_atom,hmix_ev_per_atom,hmix_mev_per_atom,hmix_ev_per_formula_unit")

for x, E_TiAlN_total in data:
    E_TiAlN_atom = E_TiAlN_total / atoms_TiAlN

    E_ref_atom = (1 - x) * E_TiN_atom + x * E_AlN_atom

    hmix_ry_atom = E_TiAlN_atom - E_ref_atom
    hmix_ev_atom = hmix_ry_atom * RY_TO_EV
    hmix_mev_atom = hmix_ev_atom * 1000
    hmix_ev_fu = hmix_ev_atom * 2

    print(f"{x},{E_TiAlN_total:.8f},{hmix_ry_atom:.10f},{hmix_ev_atom:.6f},{hmix_mev_atom:.3f},{hmix_ev_fu:.6f}")
