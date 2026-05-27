from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd

mapp = Path(__file__).resolve().parents[1]
csv_fil = mapp / "final_fixed_sqs_hmix_results.csv"
bild_fil = mapp / "final_fixed_sqs_hmix_with_reference.png"

data = pd.read_csv(csv_fil)

x = [0] + data["x"].tolist() + [1]
hmix = [0] + data["hmix_mev_per_atom"].tolist() + [0]

x_ref = [0, 0.10, 0.25044, 0.35006, 0.49950, 0.65093, 0.75156, 0.90, 1]
hmix_ref = [0, 24.502, 59.826, 81.716, 106.59, 120.52, 116.04, 75.746, 0]

plt.figure()
plt.plot(x, hmix, marker="o", label="Detta arbete")
plt.plot(x_ref, hmix_ref, marker="s", label="Referens")
plt.title("Blandningsentalpi för TiAlN")
plt.xlabel("Al-halt x")
plt.ylabel("Delta Hmix [meV/atom]")
plt.xlim(0, 1)
plt.grid(True)
plt.legend()
plt.tight_layout()
plt.savefig(bild_fil, dpi=300)

print(f"Skrev {bild_fil}")