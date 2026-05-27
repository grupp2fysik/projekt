RY_TO_EV = 13.605693122994

E_TiN_total = -848.67899733
E_ZrN_total = -1342.76289976

fu_TiN = 4
fu_ZrN = 4
fu_TiZrN = 32

data = [
    (0.25, -7777.37797040),
    (0.50, -8765.48995550),
    (0.75, -9753.74485080),
]

E_TiN_fu = E_TiN_total / fu_TiN
E_ZrN_fu = E_ZrN_total / fu_ZrN

print("x,H_mix_eV_per_atom,H_mix_meV_per_atom,H_mix_eV_per_formula_unit,H_mix_meV_per_formula_unit")
for x, E_total in data:
    E_fu = E_total / fu_TiZrN
    hmix_ev_fu = (E_fu - ((1 - x) * E_TiN_fu + x * E_ZrN_fu)) * RY_TO_EV
    hmix_ev_atom = hmix_ev_fu / 2
    print(f"{x},{hmix_ev_atom:.6f},{hmix_ev_atom*1000:.3f},{hmix_ev_fu:.6f},{hmix_ev_fu*1000:.3f}")
