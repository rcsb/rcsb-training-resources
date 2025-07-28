from rcsbapi.data import DataQuery as Query
from rcsbapi.data import ALL_STRUCTURES
from typing import Optional, List


def get_structures_with_specific_ccids(ccids: Optional[List[str]] = None):
    # Initialize query for all structures
    query = Query(
        input_type="entries",
        input_ids=ALL_STRUCTURES,
        return_data_list=[
            "rcsb_id",
            "polymer_entities.chem_comp_nstd_monomers.chem_comp.id",
            "branched_entities.chem_comp_monomers.chem_comp.id",
            "nonpolymer_entities.nonpolymer_comp.chem_comp.id"
        ]
    )

    result = query.exec(progress_bar=True)
    entries = result.get("data", {}).get("entries", [])

    chem_comp_to_pdb_map = {}

    for entry in entries:
        pdb_id = entry.get("rcsb_id")

        # Branched entities
        for branched in entry.get("branched_entities", []):
            for mono in branched.get("chem_comp_monomers", []):
                chem_id = mono.get("chem_comp", {}).get("id")
                if not chem_id:
                    continue
                if (ccids is None) or (chem_id in ccids):
                    chem_comp_to_pdb_map.setdefault(chem_id, set()).add(pdb_id)

        # Polymer entities
        for polymer in entry.get("polymer_entities", []):
            for mono in polymer.get("chem_comp_nstd_monomers", []):
                chem_id = mono.get("chem_comp", {}).get("id")
                if not chem_id:
                    continue
                if (ccids is None) or (chem_id in ccids):
                    chem_comp_to_pdb_map.setdefault(chem_id, set()).add(pdb_id)

        # Nonpolymer entities
        for nonpolymer in entry.get("nonpolymer_entities", []):
            chem_id = nonpolymer.get("nonpolymer_comp", {}).get("chem_comp", {}).get("id")
            if not chem_id:
                continue
            if (ccids is None) or (chem_id in ccids):
                chem_comp_to_pdb_map.setdefault(chem_id, set()).add(pdb_id)

    # Write to file
    output_file = "ccid-to-pdb.tsv"
    with open(output_file, "w", encoding="utf-8") as f:
        for chem_id, pdb_ids in sorted(chem_comp_to_pdb_map.items()):
            f.write(f"{chem_id}\t{' '.join(sorted(pdb_ids))}\n")

    return output_file


# Example usage:
# For specific CCIDs:
output_path = get_structures_with_specific_ccids()
# For all CCIDs:
# output_path = get_structures_with_specific_ccids()

print(f"File saved at: {output_path}")
