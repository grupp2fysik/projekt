import numpy as np
import matplotlib.pyplot as plt

# Läs DOS-fil
data = np.loadtxt("TiAlN_x05.dos", comments="#")

energy = data[:,0]
dos = data[:,1]

EF = 15.007
energy_shifted = energy - EF

plt.figure(figsize=(8,5))

plt.plot(energy_shifted, dos)
plt.axvline(0, linestyle='--')

plt.fill_between(
    energy_shifted,
    dos,
    where=(energy_shifted < 0),
    alpha=0.3
)

plt.xlabel("Energy relative to $E_F$ (eV)")
plt.ylabel("DOS")
plt.title("Full DOS for TiAlN x=0.5")

plt.grid(True)

plt.savefig("TiAlN_x05_DOS_full.png", dpi=300)

# =========================
# Zoom near EF
# =========================
plt.figure(figsize=(8,5))

plt.plot(energy_shifted, dos)
plt.axvline(0, linestyle='--')

plt.fill_between(
    energy_shifted,
    dos,
    where=(energy_shifted < 0),
    alpha=0.3
)

plt.xlabel("Energy relative to $E_F$ (eV)")
plt.ylabel("DOS")
plt.title("Zoomed DOS near $E_F$ for TiAlN x=0.5")

plt.xlim(-10,5)

plt.grid(True)

plt.savefig("TiAlN_x05_DOS_zoom.png", dpi=300)

plt.show()