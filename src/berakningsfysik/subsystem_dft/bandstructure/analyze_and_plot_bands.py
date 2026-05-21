import numpy as np
import matplotlib.pyplot as plt


def plot_bands(input_file, output_file, title, energy_zero):
    number_of_bands = 80

    data = np.loadtxt(input_file)

    points_per_band = len(data) // number_of_bands

    k = data[:points_per_band, 0]
    energies = data[:, 1].reshape(number_of_bands, points_per_band)

    labels = ["Γ", "X", "W", "K", "Γ", "L", "U", "W", "L", "K"]
    positions = [k[i] for i in [0, 20, 40, 60, 80, 100, 120, 140, 160, 180]]

    for band in energies:
        plt.plot(k, band - energy_zero, linewidth=1, color="black")

    plt.axhline(0, linestyle="--", linewidth=0.8)
    plt.xticks(positions, labels)
    plt.xlabel("Väg i Brillouin-zonen")
    plt.ylabel("Energi relativt referensnivå (eV)")
    plt.title(title)
    plt.tight_layout()
    plt.savefig(output_file, dpi=300)
    plt.close()


plot_bands(
    input_file="TiN_bands.dat.gnu",
    output_file="TiN_bandstructure.png",
    title="Kubisk TiN bandstruktur",
    energy_zero=17.0299,
)

plot_bands(
    input_file="AlN_bands.dat.gnu",
    output_file="AlN_bandstructure.png",
    title="Kubisk AlN bandstruktur",
    energy_zero=13.2165,
)

print("Created TiN_bandstructure.png")
print("Created AlN_bandstructure.png")