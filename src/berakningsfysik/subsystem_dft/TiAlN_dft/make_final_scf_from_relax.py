from pathlib import Path
import re

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

ECUTWFC = 90
ECUTRHO = 1080
KGRID = (10, 10, 10)

SYSTEMS = {
    "TiN_MP": {
        "relax_out": "convergence/relax/TiN_MP/results/vc_relax_TiN_MP.out",
        "outdir": "convergence/final_scf/TiN_MP",
        "input_name": "scf_TiN_MP_relaxed.in",
        "prefix": "TiN_MP_final_scf",
        "job_name": "TiN_scf",
    },
    "AlN_MP": {
        "relax_out": "convergence/relax/AlN_MP/results/vc_relax_AlN_MP.out",
        "outdir": "convergence/final_scf/AlN_MP",
        "input_name": "scf_AlN_MP_relaxed.in",
        "prefix": "AlN_MP_final_scf",
        "job_name": "AlN_scf",
    },
}


def extract_final_coordinates(relax_out):
    text = Path(relax_out).read_text(errors="ignore")

    match = re.search(
        r"Begin final coordinates(.*?)End final coordinates",
        text,
        flags=re.DOTALL,
    )
    if not match:
        raise RuntimeError(f"Could not find final coordinates in {relax_out}")

    block = match.group(1).strip().splitlines()

    cell_header = None
    cell_lines = []
    pos_header = None
    pos_lines = []

    i = 0
    while i < len(block):
        line = block[i].strip()

        if line.startswith("CELL_PARAMETERS"):
            cell_header = line
            cell_lines = [
                block[i + 1].strip(),
                block[i + 2].strip(),
                block[i + 3].strip(),
            ]
            i += 4
            continue

        if line.startswith("ATOMIC_POSITIONS"):
            pos_header = line
            i += 1
            while i < len(block):
                candidate = block[i].strip()
                if not candidate:
                    break
                if candidate.startswith("End"):
                    break
                # Stop if another QE section unexpectedly starts
                if candidate.startswith("CELL_PARAMETERS"):
                    break
                pos_lines.append(candidate)
                i += 1
            continue

        i += 1

    if cell_header is None or not cell_lines:
        raise RuntimeError(f"Could not parse CELL_PARAMETERS in {relax_out}")
    if pos_header is None or not pos_lines:
        raise RuntimeError(f"Could not parse ATOMIC_POSITIONS in {relax_out}")

    return cell_header, cell_lines, pos_header, pos_lines


def species_from_positions(pos_lines):
    species = []
    for line in pos_lines:
        el = line.split()[0]
        if el not in species:
            species.append(el)
    return species


def write_scf_input(system_name, data):
    relax_out = data["relax_out"]
    outdir = Path(data["outdir"])
    outdir.mkdir(parents=True, exist_ok=True)

    cell_header, cell_lines, pos_header, pos_lines = extract_final_coordinates(relax_out)
    species = species_from_positions(pos_lines)

    nat = len(pos_lines)
    ntyp = len(species)

    input_path = outdir / data["input_name"]

    with open(input_path, "w") as f:
        f.write(f"""&CONTROL
    calculation = 'scf'
    prefix = '{data["prefix"]}'
    outdir = './out'
    pseudo_dir = '{PSEUDO_DIR}'
    restart_mode = 'from_scratch'
    tstress = .true.
    tprnfor = .true.
/

&SYSTEM
    ibrav = 0
    nat = {nat}
    ntyp = {ntyp}
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

        f.write(f"\n{cell_header}\n")
        for line in cell_lines:
            f.write(line + "\n")

        f.write(f"\n{pos_header}\n")
        for line in pos_lines:
            f.write(line + "\n")

        f.write("\nK_POINTS automatic\n")
        f.write(f"{KGRID[0]} {KGRID[1]} {KGRID[2]} 0 0 0\n")

    return input_path


def write_job_script(data):
    outdir = Path(data["outdir"])
    input_name = data["input_name"]
    output_name = input_name.replace(".in", ".out")
    job_path = outdir / "job_final_scf"

    with open(job_path, "w") as f:
        f.write(f"""#!/bin/bash
#SBATCH -A liu-compute-2026-1
#SBATCH --nodes=1
#SBATCH --tasks-per-node=4
#SBATCH --time=01:00:00
#SBATCH --job-name={data["job_name"]}

module load QuantumESPRESSO/7.1-nsc1-intel-2018b-eb
export OMP_NUM_THREADS=1

mkdir -p out results

mpprun pw.x -nk 2 -in {input_name} > results/{output_name}
""")

    return job_path


def main():
    for system_name, data in SYSTEMS.items():
        input_path = write_scf_input(system_name, data)
        job_path = write_job_script(data)

        print("=" * 60)
        print(system_name)
        print(f"Wrote input: {input_path}")
        print(f"Wrote job:   {job_path}")
        print("Run:")
        print(f"cd {data['outdir']}")
        print("sbatch job_final_scf")


if __name__ == "__main__":
    main()