from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd

mapp = Path(__file__).resolve().parent

egen_csv = mapp / "TiZrN_hmix_results.csv"
referens_csv = mapp / "liu_tizrn_reference_digitized.csv"
bild_fil = mapp / "TiZrN_hmix_with_reference.png"

egen_data = pd.read_csv(egen_csv)
referens_data = pd.read_csv(referens_csv)

x = [0] + egen_data["x"].tolist() + [1]
hmix = [0] + egen_data["H_mix_meV_per_formula_unit"].tolist() + [0]

plt.figure(figsize=(5, 5))
plt.plot(x, hmix, marker="o", label="Detta arbete")
plt.plot(
    referens_data["x"],
    referens_data["hmix_mev_per_formula_unit"],
    marker="s",
    label="Referens",
)

plt.title("Blandningsentalpi för TiZrN")
plt.xlabel("Zr-halt x")
plt.ylabel("Delta Hmix [meV/formelenhet]")
plt.xlim(0, 1)
plt.ylim(-20, 160)
plt.xticks([0, 0.25, 0.5, 0.75, 1], ["TiN", "0.25", "0.5", "0.75", "ZrN"])
plt.grid(True)
plt.legend()
plt.tight_layout()
plt.savefig(bild_fil, dpi=300)

print(f"Skrev {bild_fil}")
