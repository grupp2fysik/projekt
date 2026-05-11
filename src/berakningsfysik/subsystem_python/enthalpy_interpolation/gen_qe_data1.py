import argparse
import numpy as np
from pathlib import Path
import textwrap

RY_TO_EV = 13.60593

DEFAULT_XS = [i / 16 for i in range(17)]

E_TIN = -3.80
E_ALN = -3.30

A_TIN = 4.24
A_ALN = 4.05

NATOMS = 8


def hmix_single_peak(x: float) -> float:
    """
    En vanlig enkelpucklig blandningsentalpi.
    Max ungefär kring x = 0.5.
    """
    z = 2.0 * x - 1.0
    return x * (1.0 - x) * (0.80 + 0.20 * z - 0.10 * z**2)


def hmix_double_peak(x: float) -> float:
    """
    Dubbelpucklig blandningsentalpi.

    Formen är fortfarande Redlich-Kister:
        H = x(1-x) * (L0 + L2 z^2)

    Med L2 > L0 får man två lokala maxima,
    symmetriskt placerade kring x = 0.5.
    """
    z = 2.0 * x - 1.0

    L0 = 0.18
    L2 = 1.00

    return x * (1.0 - x) * (L0 + L2 * z**2)


def hmix_eva(x: float, peaks: int) -> float:
    """
    Välj syntetisk blandningsentalpi.
    """
    if peaks == 1:
        return hmix_single_peak(x)

    if peaks == 2:
        return hmix_double_peak(x)

    raise ValueError(f"Stödjer bara peaks=1 eller peaks=2, fick {peaks}")


def lattice_a(x: float) -> float:
    """
    Syntetisk gitterparameter i Å.

    Viktigt: detta ska vara en längd, inte en energi.
    """
    bowing = 0.03 * x * (1.0 - x)
    return (1.0 - x) * A_TIN + x * A_ALN + bowing


def energy_per_atom_eV(x: float, peaks: int) -> float:
    return (1.0 - x) * E_TIN + x * E_ALN + hmix_eva(x, peaks)


def total_energy_ry(x: float, peaks: int) -> float:
    return energy_per_atom_eV(x, peaks) * NATOMS / RY_TO_EV


def composition_counts(x: float, metal_sites: int = 4) -> tuple[int, int, int]:
    n_al = int(round(x * metal_sites))
    n_ti = metal_sites - n_al
    n_n = 4
    return n_ti, n_al, n_n


def make_qe_out_text(x: float, peaks: int) -> str:
    etot_ry = total_energy_ry(x, peaks)
    a = lattice_a(x)
    n_ti, n_al, n_n = composition_counts(x)

    a1 = (a, 0.0, 0.0)
    a2 = (0.0, a * (1.0 + 0.002), 0.0)
    a3 = (0.0, 0.0, a * (1.0 - 0.0015))

    return textwrap.dedent(
        f"""\
        Program PWSCF v.7.2 starts on  01-Apr-2026 at 10:00:00

             number of atoms/cell      =            {NATOMS}
             number of atomic types    =            3
             number of electrons       =         64.00
             lattice parameter (alat)  =       8.00000000  a.u.
             celldm(1)=   8.00000000

        bravais-lattice index     =            0
        unit-cell volume          =      512.0000 (a.u.)^3
        kinetic-energy cutoff     =       60.0000  Ry
        charge density cutoff     =      480.0000  Ry

        atomic species   valence    mass     pseudopotential
        Ti               12.00    47.86700   Ti.UPF
        Al                3.00    26.98154   Al.UPF
        N                 5.00    14.00670   N.UPF

        Self-consistent Calculation

             iteration #  1     ecut=    60.00 Ry     beta= 0.70
             total cpu time spent up to now is        2.0 secs
             estimated scf accuracy    <       0.01000000 Ry
             total energy              =     {etot_ry + 0.05000000: .8f} Ry

             iteration #  2     ecut=    60.00 Ry     beta= 0.70
             total cpu time spent up to now is        3.0 secs
             estimated scf accuracy    <       0.00100000 Ry
             total energy              =     {etot_ry + 0.00500000: .8f} Ry

             iteration #  3     ecut=    60.00 Ry     beta= 0.70
             total cpu time spent up to now is        4.0 secs
             estimated scf accuracy    <       0.00001000 Ry
             total energy              =     {etot_ry + 0.00050000: .8f} Ry

        convergence has been achieved in   9 iterations

        End of self-consistent calculation

        k = 0.0000 0.0000 0.0000 ( 14712 PWs)   bands (ev):

        -17.3307  -9.3182  -9.3176  -9.3173

        !    total energy              =     {etot_ry: .8f} Ry
             estimated scf accuracy    <       0.00000064 Ry

        CELL_PARAMETERS (angstrom)
          {a1[0]: .8f}  {a1[1]: .8f}  {a1[2]: .8f}
          {a2[0]: .8f}  {a2[1]: .8f}  {a2[2]: .8f}
          {a3[0]: .8f}  {a3[1]: .8f}  {a3[2]: .8f}

        ATOMIC_POSITIONS (angstrom)
        Ti  0.000000  0.000000  0.000000
        Ti  0.500000  0.500000  0.000000
        Al  0.500000  0.000000  0.500000
        N   0.000000  0.500000  0.500000
        N   0.250000  0.250000  0.250000
        N   0.750000  0.250000  0.250000
        N   0.250000  0.750000  0.250000
        N   0.250000  0.250000  0.750000

        Synthetic composition summary:
          x = {x:.3f}
          Ti atoms = {n_ti}
          Al atoms = {n_al}
          N atoms = {n_n}
          H_mix = {hmix_eva(x, peaks):.8f} eV/atom
        """
    )


def parse_x_grid(text: str) -> list[float]:
    return [float(v.strip()) for v in text.split(",") if v.strip()]


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Skapa syntetiska Quantum ESPRESSO .out-filer."
    )

    parser.add_argument(
        "--out-dir",
        default="qe_outputs/qe_outputs2",
        help="Mapp där .out-filerna sparas.",
    )

    parser.add_argument(
        "--peaks",
        type=int,
        choices=[1, 2],
        default=1,
        help="Antal lokala maxima i syntetisk blandningsentalpi. Välj 1 eller 2.",
    )

    parser.add_argument(
        "--xs",
        default=",".join(f"{x:.4f}" for x in DEFAULT_XS),
        help="Kommaseparerad lista med x-värden.",
    )

    args = parser.parse_args()

    xs = parse_x_grid(args.xs)

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    for x in xs:
        filename = out_dir / f"TiAlN_x={x:.4f}.out"
        filename.write_text(make_qe_out_text(x, args.peaks), encoding="utf-8")

    print(f"Skapade {len(xs)} syntetiska .out-filer i: {out_dir.resolve()}")
    print(f"Vald form: peaks={args.peaks}")

    if args.peaks == 2:
        print("Dubbelpucklig ΔH_mix aktiverad.")
        print("Förväntade maxima ligger ungefär nära x ≈ 0.18 och x ≈ 0.82.")


if __name__ == "__main__":
    main()