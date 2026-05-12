from ase.build import bulk
from icet import ClusterSpace
from icet.tools.structure_generation import (generate_sqs,
                                             generate_sqs_from_supercells,
                                             generate_sqs_by_enumeration,
                                             generate_target_structure)
from ase.io import write
from ase.io import read
from icet.tools.structure_generation import compare_cluster_vectors

def generate_sqs_supercell(x_al=0.5, supercell_size= (4,4,4),cubicfalsetrue=True, filename = "sqs_supercell.cif"):
    #Uträkning av a:
    a= 4.07*x_al+(1-x_al)*4.24
    # 1. Struktur
    structure = bulk('TiN', crystalstructure='rocksalt', a=a,cubic= cubicfalsetrue) #säger strukturen hos det vi ska skapa. I detta fall är det rocksalt och är TiN-struktur.
    # 2. Sublattices
    if(cubicfalsetrue):
        chemical_symbols = [
            ['Ti', 'Al'],  # atom 0 är Ti-plats
            ['N'],         # atom 1 är N-plats
            ['Ti', 'Al'],  # atom 2 är Ti-plats
            ['N'],         # atom 3 är N-plats
            ['Ti', 'Al'],  # atom 4 är Ti-plats
            ['N'],         # atom 5 är N-plats
            ['Ti', 'Al'],  # atom 6 är Ti-plats
            ['N']          # atom 7 är N-plats
        ]
    else:
        chemical_symbols = [['Ti', 'Al'], ['N']] #säger att det är en sublattice och att det är titan och aluminium i ena och N i andra.
    

    # 3. Cluster space
    cutoff = 1.5* a #Hur många grannskal kommer beräknas (i ångström)
    cs = ClusterSpace(structure, [cutoff], chemical_symbols=chemical_symbols) #första är strukturen, andra är hur många grannskal (desto högre desto bättre men mer beräkningstid), 
    # Tredje är sublatticesen och vad som är i dem.

    # 4. Target concentrations, säger procenten av sublattices de ska vara
    target_concentrations = {
        'A': {'Ti': 1 - x_al, 'Al': x_al},
        'B': {'N': 1.0}
    }

    supercells = [structure.repeat(supercell_size)]
    sqs_supercell = generate_sqs_from_supercells(cluster_space=cs,
                                   supercells=supercells,
                                   n_steps=160000,
                                   target_concentrations=target_concentrations)

    cv_sqs = cs.get_cluster_vector(sqs_supercell)

    print("Cluster vector för SQS:")
    print(cv_sqs)
    
    write(filename, sqs_supercell)

    return sqs_supercell


def cif_to_qe_input(cif_file="test_supercell.cif",  # Namnet på inputfilen (.cif med struktur)
                    output_file="qe_input.in",      # Namnet på outputfilen (QE input)
                    pseudo_dir="./pseudo/"):        # Sökväg till pseudopotentialer (används av QE senare)

    # Läser in strukturen från .cif-filen och skapar ett ASE Atoms-objekt
    atoms = read(cif_file)

    symbols = atoms.get_chemical_symbols()
    unique_symbols = sorted(set(symbols))

    # Öppnar (eller skapar) outputfilen och skriver till den
    with open(output_file, "w") as f:

        

        # ATOMIC SPECIES
        # Kopplar varje element till massa + pseudopotential
        pseudo_map = {
            "Ti": ("47.867", "Ti.pbe-spn-kjpaw_psl.1.0.0.UPF"),
            "Al": ("26.982", "Al.pbe-n-kjpaw_psl.1.0.0.UPF"),
            "N":  ("14.007", "N.pbe-n-kjpaw_psl.1.0.0.UPF"),
        }

        f.write("ATOMIC_SPECIES\n")
        for s in unique_symbols:
            mass, pseudo = pseudo_map[s]
            # Skriver element, atomvikt och pseudopotential-fil
            f.write(f"{s} {mass} {pseudo}\n")

        f.write("\n")

        # CELL PARAMETERS
        # Skriver ut enhetscellens vektorer (i angstrom)
        f.write("CELL_PARAMETERS angstrom\n")
        for v in atoms.get_cell():
            f.write(f"{v[0]} {v[1]} {v[2]}\n")

        f.write("\n")

        # ATOMIC POSITIONS
        # Skriver alla atomers positioner i cellen
        f.write("ATOMIC_POSITIONS angstrom\n")
        for s, pos in zip(symbols, atoms.get_positions()):
            f.write(f"{s} {pos[0]} {pos[1]} {pos[2]}\n")

        f.write("\n")

        # K-POINTS
        f.write("K_POINTS automatic\n")
        f.write("2 2 2 0 0 0\n")

    # Bekräftelse att filen skapats
    print(f" input skapad: {output_file}")

if __name__ == "__main__":
    x_al = 0.875
    label = "x0875"

    cif_file = f"TiAlN_{label}_sqs_64atoms.cif"
    qe_file = f"TiAlN_{label}_sqs_64atoms_qe_structure.in"

    sc = generate_sqs_supercell(
        x_al=x_al,
        supercell_size=(2, 2, 2),
        cubicfalsetrue=True,
        filename=cif_file
    )

    print("Klar!")
    print("Antal atomer:", len(sc))

    cif_to_qe_input(
        cif_file=cif_file,
        output_file=qe_file
    )