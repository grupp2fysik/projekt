import os
from mp_api.client import MPRester
from pymatgen.io.cif import CifWriter

API_KEY = os.environ.["MP_API_KEY"]

materials = {
    "TiN_mp-492": "mp-492",
    "AlN_mp-1330": "mp-1330"
}

with MPRester(API_KEY) as mpr:
    for name, mpid in materials.items():
        docs = mpr.materials.summary.search(
            material_ids=[mpid],
            fields=["structure", "formula_pretty", "symmetry"]
        )

        if len(docs) != 1:
            raise RuntimeError(f"Expected one result for {mpid}, got {len(docs)}")

        structure = docs[0].structure

        # Gör en konventionell cell för tydlig rapportering och enklare QE-input.
        conv = structure.to_conventional()

        outname = f"structures/{name}_conventional.cif"
        CifWriter(conv).write_file(outname)

        print("=" * 60)
        print(name)
        print("MP-ID:", mpid)
        print("formula:", docs[0].formula_pretty)
        print("symmetry:", docs[0].symmetry)
        print("sites in conventional cell:", len(conv))
        print("lattice:")
        print(conv.lattice)
        print("wrote:", outname)