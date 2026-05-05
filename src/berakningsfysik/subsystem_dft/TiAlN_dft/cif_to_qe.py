from pathlib import Path
from pymatgen.core import Structure

# ================================
# User settings
# ================================

PSEUDO_DIR = "../../../pseudos/pseudos_pslib"

PSEUDOS = {
    "Ti": "Ti.pbe-spn-kjpaw_psl.1.0.0.UPF",
    "Al": "Al.pbe-n-kjpaw_psl.1.0.0.UPF",
    "N":  "N.pbe-n-kjpaw_psl.1.0.0.UPF",
}

MASSES = {
    "Ti": 47.867,
    "Al": 26.982,
    "N": 14.007,
}

# Good starting values for PSLibrary PAW.
# These are not final; they will be convergence tested.
ECUTWFC = 80
ECUTRHO = 800

KGRID = (10, 10, 10)

STRUCTURES = {
    "TiN": {
        "cif": "structures/TiN_mp-492_conventional.cif",
        "outdir": "convergence/cutoff/TiN_MP_test",
        "prefix": "TiN_MP_test",
        "input_name": "scf_TiN_MP_test.in",
    },
    "AlN": {
        "cif": "structures/AlN_mp-1330_conventional.cif",
        "outdir": "convergence/cutoff/AlN_MP_test",
        "prefix": "AlN_MP_test",
        "input_name": "scf_AlN_MP_test.in",
    },
}

SPECIES_ORDER = ["Ti", "Al", "N"]


def unique_species_in_order(structure):
    present = {str(site.specie) for site in structure}
    return [el for el in SPECIES_ORDER if el in present]


def write_qe_input(cif_path, out_path, prefix):
    structure = Structure.from_file(cif_path)
    species = unique_species_in_order(structure)

    missing_pseudos = [el for el in species if el not in PSEUDOS]
    if missing_pseudos:
        raise ValueError(f"Missing pseudopotentials for: {missing_pseudos}")

    out_path.parent.mkdir(parents=True, exist_ok=True)

    with open(out_path, "w") as f:
        f.write(f"""&CONTROL
    calculation = 'scf'
    prefix = '{prefix}'
    outdir = './out'
    pseudo_dir = '{PSEUDO_DIR}'
    restart_mode = 'from_scratch'
    tstress = .true.
    tprnfor = .true.
/

&SYSTEM
    ibrav = 0
    nat = {len(structure)}
    ntyp = {len(species)}
    ecutwfc = {ECUTWFC}
    ecutrho = {ECUTRHO}
    occupations = 'smearing'
    smearing = 'gaussian'
    degauss = 0.005
/

&ELECTRONS
    conv_thr = 1.0d-8
    mixing_beta = 0.2
    diagonalization = 'david'
/

ATOMIC_SPECIES
""")

        for el in species:
            f.write(f"{el:<2} {MASSES[el]:>8.3f} {PSEUDOS[el]}\n")

        f.write("\nCELL_PARAMETERS angstrom\n")
        for vec in structure.lattice.matrix:
            f.write(f"{vec[0]: .12f} {vec[1]: .12f} {vec[2]: .12f}\n")

        f.write("\nATOMIC_POSITIONS crystal\n")
        for site in structure:
            x, y, z = site.frac_coords
            # Wrap coordinates into [0, 1) to avoid tiny negative numerical noise.
            x, y, z = x % 1.0, y % 1.0, z % 1.0
            f.write(f"{str(site.specie):<2} {x: .12f} {y: .12f} {z: .12f}\n")

        f.write("\nK_POINTS automatic\n")
        f.write(f"{KGRID[0]} {KGRID[1]} {KGRID[2]} 0 0 0\n")


def write_job_script(job_path, input_name, output_name, job_name):
    with open(job_path, "w") as f:
        f.write(f"""#!/bin/bash
#SBATCH -A liu-compute-2026-1
#SBATCH --nodes=1
#SBATCH --tasks-per-node=4
#SBATCH --time=00:30:00
#SBATCH --job-name={job_name}

module load QuantumESPRESSO/7.1-nsc1-intel-2018b-eb
export OMP_NUM_THREADS=1

mkdir -p out results

mpprun pw.x -nk 2 -in {input_name} > results/{output_name}
""")


def main():
    for material, data in STRUCTURES.items():
        cif_path = Path(data["cif"])
        outdir = Path(data["outdir"])
        input_path = outdir / data["input_name"]
        job_path = outdir / "job_qe"
        output_name = data["input_name"].replace(".in", ".out")

        if not cif_path.exists():
            raise FileNotFoundError(f"Could not find CIF file: {cif_path}")

        write_qe_input(
            cif_path=cif_path,
            out_path=input_path,
            prefix=data["prefix"],
        )

        write_job_script(
            job_path=job_path,
            input_name=data["input_name"],
            output_name=output_name,
            job_name=f"{material}_MP_test",
        )

        print("=" * 60)
        print(f"Wrote QE input for {material}")
        print(f"Input: {input_path}")
        print(f"Job:   {job_path}")
        print(f"Run on Sigma with:")
        print(f"cd {outdir}")
        print("sbatch job_qe")


if __name__ == "__main__":
    main()