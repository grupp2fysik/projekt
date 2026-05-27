import pandas as pd
import matplotlib.pyplot as plt

df = pd.read_csv("TiZrN_hmix_results.csv")

x = [0] + df["x"].tolist() + [1]
h = [0] + df["H_mix_meV_per_formula_unit"].tolist() + [0]

plt.figure(figsize=(5, 5))
plt.plot(x, h, marker="o")
plt.title("Blandningsentalpi för TiZrN")
plt.xlabel("Zr-halt x")
plt.ylabel("Blandningsentalpi [meV/formelenhet]")
plt.xlim(0, 1)
plt.ylim(0, 140)
plt.xticks([0, 0.25, 0.5, 0.75, 1], ["TiN", "0.25", "0.5", "0.75", "ZrN"])
plt.tight_layout()
plt.savefig("TiZrN_hmix_mev_fu.png", dpi=300)
plt.show()