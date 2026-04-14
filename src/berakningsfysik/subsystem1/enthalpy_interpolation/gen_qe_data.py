from __future__ import annotations

from pathlib import Path
import textwrap

RY_TO_EV = 13.60593

# Molbråk för Al 
xs = [0.0, 0.125, 0.25, 0.375, 0.5, 0.625, 0.75, 0.75, 0.875, 1.0]

# Syntetiska referensenergier i eV/atom
E_TIN = -3.80
E_ALN = -3.30

NATOMS = 8

def hmix_eva(x: float) -> float:
    """Syntetisk blandningsentalpi i eV/atom."""
    z = 2.0 * x - 1.0
    return x * (1.0 - x) * (0.80 + 0.20 * z - 0.10 * z**2)


def lattice_a(x: float) -> float:
    """Syntetisk gitterparameter i Ångström."""
    return (1.0 - x) * E_TIN + x * E_ALN + hmix_eva(x)


def energy_per_atom_eV(x: float) -> float:
    return (1.0 - x) * E_TIN + x * E_ALN + hmix_eva(x)


def total_energy_ry(x: float) -> float:
    return energy_per_atom_eV(x) * NATOMS / RY_TO_EV


def composition_counts(x: float, metal_sites: int = 4) -> tuple[int, int, int]:
    n_al = int(round(x * metal_sites))
    n_ti = metal_sites - n_al
    n_n = 4
    return n_ti, n_al, n_n


def make_qe_out_text(x: float) -> str:
    etot_ry = total_energy_ry(x)
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
        """
    )


def main() -> None:
    out_dir = Path("qe_outputs")
    out_dir.mkdir(exist_ok=True)

    for x in xs:
        filename = out_dir / f"TiAlN_x={x:.3f}.out"
        filename.write_text(make_qe_out_text(x), encoding="utf-8")

    print(f"Skapade {len(xs)} syntetiska .out-filer i: {out_dir.resolve()}")


if __name__ == "__main__":
    main()