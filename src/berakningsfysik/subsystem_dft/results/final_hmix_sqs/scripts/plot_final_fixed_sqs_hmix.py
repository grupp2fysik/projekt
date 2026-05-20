from pathlib import Path

import pandas as pd
import matplotlib.pyplot as plt

mapp = Path(__file__).resolve().parents[1]
csv_fil = mapp / "final_fixed_sqs_hmix_results.csv"
bild_fil = mapp / "final_fixed_sqs_hmix_plot.png"

data = pd.read_csv(csv_fil)

plt.figure()
plt.plot(data["x"], data["hmix_ev_per_atom"], marker="o")
plt.title("Blandningsentalpi för TiAlN")
plt.xlabel("Al-halt x")
plt.ylabel("Delta Hmix [eV/atom]")
plt.xlim(0.0, 1.0)
plt.ylim(-0.10, 0.25)
plt.grid(True)
plt.tight_layout()
plt.savefig(bild_fil, dpi=300)

print(f"Skrev {bild_fil}")
