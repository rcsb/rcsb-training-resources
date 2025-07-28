from rcsbapi.data import DataQuery as Query
from rcsbapi.data import ALL_STRUCTURES
from typing import Optional, List


def get_entries_by_ccids_with_entity_check(ccids: Optional[List[str]] = None):
    """
    Retrieves PDB entries that contain specific chemical component IDs (ccids),
    or all chemical components if ccids is None. Uses RCSB Search API with strict
    nested entity filtering (polymer, branched, nonpolymer).
    """

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
        branched_entities = entry.get("branched_entities")
        if branched_entities:
            for branched in branched_entities:
                chem_comp_monomers = branched.get("chem_comp_monomers")
                if chem_comp_monomers:
                    for mono in chem_comp_monomers:
                        chem_comp = mono.get("chem_comp")
                        if chem_comp:
                            chem_id = chem_comp.get("id")
                            if ccids is None or chem_id in ccids:
                                chem_comp_to_pdb_map.setdefault(chem_id, set()).add(pdb_id)

        # Polymer entities
        polymer_entities = entry.get("polymer_entities")
        if polymer_entities:
            for polymer in polymer_entities:
                chem_comp_nstd_monomers = polymer.get("chem_comp_nstd_monomers")
                if chem_comp_nstd_monomers:
                    for mono in chem_comp_nstd_monomers:
                        chem_comp = mono.get("chem_comp")
                        if chem_comp:
                            chem_id = chem_comp.get("id")
                            if ccids is None or chem_id in ccids:
                                chem_comp_to_pdb_map.setdefault(chem_id, set()).add(pdb_id)

        # Nonpolymer entities
        nonpolymer_entities = entry.get("nonpolymer_entities")
        if nonpolymer_entities:
            for nonpolymer in nonpolymer_entities:
                nonpolymer_comp = nonpolymer.get("nonpolymer_comp")
                if nonpolymer_comp:
                    chem_comp = nonpolymer_comp.get("chem_comp")
                    if chem_comp:
                        chem_id = chem_comp.get("id")
                        if ccids is None or chem_id in ccids:
                            chem_comp_to_pdb_map.setdefault(chem_id, set()).add(pdb_id)

    # Write results to file
    output_file = "ccid-to-pdb.tsv"
    with open(output_file, "w", encoding="utf-8") as f:
        for chem_id, pdb_ids in sorted(chem_comp_to_pdb_map.items()):
            f.write(f"{chem_id}\t{' '.join(sorted(pdb_ids))}\n")

    return output_file


# Example usage:

# Run with specific CCIDs
# output_path = get_entries_by_ccids_with_entity_check(["ATP", "GTP", "FAD"])

# Run for all CCIDs in all entries
output_path = get_entries_by_ccids_with_entity_check()

print(f"File saved at: {output_path}")
