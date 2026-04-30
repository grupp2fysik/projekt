from pathlib import Path
from pymatgen.core import Structure

for f in Path("structures").glob("*.cif"):
    s = Structure.from_file(f)
    print("=" * 60)
    print(f)
    print("formula:", s.composition.reduced_formula)
    print("sites:", len(s))
    print("lattice:")
    print(s.lattice)
    print("species:", sorted(set(str(site.specie) for site in s)))
