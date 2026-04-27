from ase.build import bulk
from icet import ClusterSpace
from icet.tools.structure_generation import (generate_sqs,
                                             generate_sqs_from_supercells,
                                             generate_sqs_by_enumeration,
                                             generate_target_structure)
from ase.io import write
from ase.io import read


def generate_tialn_sqs(a=4.24, x_al=0.5, size=8, filename="sqs.cif"): #Första är vet inte exakt, andra är hur många atomer på enetscellen, tredje är namnet på filen som skapas
#INPUTS: x_al= proportionen av Al i legeringen, size = antal atomer, filename=namn på filen som skapas.

    # 1. Struktur
    structure = bulk('TiN', crystalstructure='rocksalt', a=a) #säger strukturen hos det vi ska skapa. I detta fall är det rocksalt och är TiN-struktur.

    # 2. Sublattices
    chemical_symbols = [['Ti', 'Al'], ['N']] #säger att det är en sublattice och att det är titan och aluminium i ena och N i andra.

    # 3. Cluster space
    cutoff = 1.5* a #Hur många grannskal kommer beräknas (i ångström)
    cs = ClusterSpace(structure, [7.0], chemical_symbols=chemical_symbols) #första är strukturen, andra är hur många grannskal (desto högre desto bättre men mer beräkningstid), 
    # Tredje är sublatticesen och vad som är i dem.

    # 4. Target concentrations, säger procenten av sublattices de ska vara
    target_concentrations = {
        'A': {'Ti': 1 - x_al, 'Al': x_al},
        'B': {'N': 1.0}
    }

    # 5. Generera SQS
    sqs = generate_sqs(
        cluster_space=cs,
        max_size=size,
        include_smaller_cells=False, #Ifall denna är true, så kollar den också 1 till max_size atomer, och jämför resultatet och tar den bästa
        target_concentrations=target_concentrations, #Vad målet av de olika sublattices koncentrationerna kommer vara.
        n_steps=10000 #Antal steg den gör beräkningar, desto mer desto bättre beräkning men mer tidskrävande. Detta kan kanske också vara en input senare
    )

    # 6. Spara
    write(filename, sqs) #Skriver in i filnamnet

    return sqs

"""def generate_sqs_supercell(a=4.24, x_al=0.5, sqs_size=8, supercell_size=(4,4,4), filename="sqs_supercell.cif"): #Den här är gammal kod, men bevarar den in case of.
    
    #Genererar först en SQS och gör sedan en supercell av den


    # 1. Skapa SQS med din gamla funktion
    sqs = generate_tialn_sqs(a=a, x_al=x_al, size=sqs_size)

    # 2. Skapa supercell
    supercell = sqs.repeat(supercell_size)

    # 3. Spara
    write(filename, supercell)

    return supercell"""

def generate_sqs_supercell(a=4.24, x_al=0.5, supercell_size= (4,4,4), filename = "sqs_supercell.cif"):
    # 1. Struktur
    structure = bulk('TiN', crystalstructure='rocksalt', a=a) #säger strukturen hos det vi ska skapa. I detta fall är det rocksalt och är TiN-struktur.

    # 2. Sublattices
    chemical_symbols = [['Ti', 'Al'], ['N']] #säger att det är en sublattice och att det är titan och aluminium i ena och N i andra.

    # 3. Cluster space
    cutoff = 1.5* a #Hur många grannskal kommer beräknas (i ångström)
    cs = ClusterSpace(structure, [7.0], chemical_symbols=chemical_symbols) #första är strukturen, andra är hur många grannskal (desto högre desto bättre men mer beräkningstid), 
    # Tredje är sublatticesen och vad som är i dem.

    # 4. Target concentrations, säger procenten av sublattices de ska vara
    target_concentrations = {
        'A': {'Ti': 1 - x_al, 'Al': x_al},
        'B': {'N': 1.0}
    }

    supercells = [structure.repeat(supercell_size)]
    sqs_supercell = generate_sqs_from_supercells(cluster_space=cs,
                                   supercells=supercells,
                                   n_steps=10000,
                                   target_concentrations=target_concentrations)
    write(filename, sqs_supercell)

    return sqs_supercell


def cif_to_qe_input(cif_file="test_supercell.cif",  # Namnet på inputfilen (.cif med struktur)
                    output_file="qe_input.in",      # Namnet på outputfilen (QE input)
                    pseudo_dir="./pseudo/"):        # Sökväg till pseudopotentialer (används av QE senare)

    # Läser in strukturen från .cif-filen och skapar ett ASE Atoms-objekt
    atoms = read(cif_file)

    # Öppnar (eller skapar) outputfilen och skriver till den
    with open(output_file, "w") as f:

        # CONTROL-block
        # Styr själva beräkningen i Quantum ESPRESSO
        f.write("&CONTROL\n")
        f.write("  calculation = 'scf'\n")  # Typ av beräkning (self-consistent field)
        f.write("  prefix = 'tialn'\n")     # Namn på beräkningen (används för outputfiler)
        f.write("  outdir = './outdir/'\n") # Var temporära filer sparas
        f.write(f"  pseudo_dir = '{pseudo_dir}'\n")  # Var pseudopotentialerna finns
        f.write("/\n\n")

        # SYSTEM-block
        # Beskriver systemets fysik och struktur
        symbols = atoms.get_chemical_symbols()      # Lista på alla atomtyper (t.ex. ['Ti', 'N', ...])
        unique_symbols = sorted(set(symbols))       # Unika atomtyper (t.ex. ['Al', 'N', 'Ti'])

        f.write("&SYSTEM\n")
        f.write("  ibrav = 0\n")                    # Vi anger cellen manuellt (CELL PARAMETERS)
        f.write(f"  nat = {len(atoms)}\n")          # Antal atomer i systemet
        f.write(f"  ntyp = {len(unique_symbols)}\n")# Antal olika atomtyper
        f.write("  ecutwfc = 200\n")                # Cutoff för vågfunktioner (viktig för noggrannhet)
        f.write("  ecutrho = 800\n")                # Cutoff för laddningstäthet
        f.write("  occupations = 'smearing'\n")     # Hantering av elektroner (bra för metaller)
        f.write("  smearing = 'gaussian'\n")        # Typ av smearing
        f.write("  degauss = 0.005\n")              # Bredd på smearing
        f.write("/\n\n")

        # ----------------------
        # ELECTRONS-block
        # ----------------------
        # Inställningar för den numeriska lösningen
        f.write("&ELECTRONS\n")
        f.write("  conv_thr = 1.0d-8\n")            # Konvergenskrav (hur noggrant lösningen ska vara)
        f.write("  mixing_beta = 0.2\n")            # Hur snabbt lösningen uppdateras
        f.write("/\n\n")

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
        f.write("4 4 4 0 0 0\n")

    # Bekräftelse att filen skapats
    print(f" input skapad: {output_file}")

if __name__ == "__main__":
    sc = generate_sqs_supercell(
        a=4.24,
        x_al=0.5,
        supercell_size=(2,2,2),
        filename="test_supercell.cif"
    )
    #sc = generate_tialn_sqs()
    print("Klar!")
    print("Antal atomer:", len(sc))
    cif_to_qe_input("test_supercell.cif")
