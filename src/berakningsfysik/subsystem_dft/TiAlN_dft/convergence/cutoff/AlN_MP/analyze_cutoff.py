from pathlib import Path
import re
import csv

RY_TO_EV = 13.605693122994
NAT = 8

rows = []

for path in Path("results").glob("scf_AlN_ecut*.out"):
    m_ecut = re.search(r"ecut(\d+)", path.name)
    if not m_ecut:
        continue

    ecut = int(m_ecut.group(1))
    text = path.read_text(errors="ignore")

    energies = re.findall(r"!\s+total energy\s+=\s+(-?\d+\.\d+)", text)
    job_done = "JOB DONE" in text

    if not energies:
        rows.append({
            "ecutwfc_Ry": ecut,
            "E_Ry": None,
            "E_eV_per_atom": None,
            "job_done": job_done,
        })
        continue

    E_Ry = float(energies[-1])
    E_eV_per_atom = E_Ry * RY_TO_EV / NAT

    rows.append({
        "ecutwfc_Ry": ecut,
        "E_Ry": E_Ry,
        "E_eV_per_atom": E_eV_per_atom,
        "job_done": job_done,
    })

rows.sort(key=lambda r: r["ecutwfc_Ry"])

valid = [r for r in rows if r["E_eV_per_atom"] is not None]
if not valid:
    raise RuntimeError("No valid energies found.")

E_ref = valid[-1]["E_eV_per_atom"]

for r in rows:
    if r["E_eV_per_atom"] is None:
        r["delta_meV_per_atom_vs_highest"] = None
    else:
        r["delta_meV_per_atom_vs_highest"] = abs(r["E_eV_per_atom"] - E_ref) * 1000

out_csv = "AlN_cutoff_convergence.csv"
with open(out_csv, "w", newline="") as f:
    writer = csv.DictWriter(
        f,
        fieldnames=[
            "ecutwfc_Ry",
            "E_Ry",
            "E_eV_per_atom",
            "delta_meV_per_atom_vs_highest",
            "job_done",
        ],
    )
    writer.writeheader()
    writer.writerows(rows)

print(f"Wrote {out_csv}")
print()
print(f"{'ecut':>8} {'E_Ry':>18} {'E_eV/atom':>18} {'delta meV/atom':>18} {'JOB DONE':>10}")
for r in rows:
    print(
        f"{r['ecutwfc_Ry']:8} "
        f"{r['E_Ry'] if r['E_Ry'] is not None else 'NA':>18} "
        f"{r['E_eV_per_atom'] if r['E_eV_per_atom'] is not None else 'NA':>18} "
        f"{r['delta_meV_per_atom_vs_highest'] if r['delta_meV_per_atom_vs_highest'] is not None else 'NA':>18} "
        f"{str(r['job_done']):>10}"
    )