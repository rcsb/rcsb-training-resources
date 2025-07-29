from rcsbapi.data import DataQuery as Query
from rcsbapi.data import ALL_STRUCTURES


def get_all_chem_comp_ids_and_write_to_file():
    # Initialize the query
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

    # Execute the query
    result = query.exec(progress_bar=True)
    entries = result.get("data", {}).get("entries", [])

    chem_comp_to_pdb_map = {}

    for entry in entries:
        pdb_id = entry.get("rcsb_id")

        # 1. Branched Entities
        branched_entities = entry.get("branched_entities")
        if branched_entities:
            for branched in branched_entities:
                chem_comp_monomers = branched.get("chem_comp_monomers")
                if chem_comp_monomers:
                    for mono in chem_comp_monomers:
                        chem_comp = mono.get("chem_comp")
                        if chem_comp:
                            chem_id = chem_comp.get("id")
                            if chem_id:
                                chem_comp_to_pdb_map.setdefault(chem_id, set()).add(pdb_id)

        # 2. Polymer Entities
        polymer_entities = entry.get("polymer_entities")
        if polymer_entities:
            for polymer in polymer_entities:
                chem_comp_nstd_monomers = polymer.get("chem_comp_nstd_monomers")
                if chem_comp_nstd_monomers:
                    for mono in chem_comp_nstd_monomers:
                        chem_comp = mono.get("chem_comp")
                        if chem_comp:
                            chem_id = chem_comp.get("id")
                            if chem_id:
                                chem_comp_to_pdb_map.setdefault(chem_id, set()).add(pdb_id)

        # 3. Nonpolymer Entities
        nonpolymer_entities = entry.get("nonpolymer_entities")
        if nonpolymer_entities:
            for nonpolymer in nonpolymer_entities:
                nonpolymer_comp = nonpolymer.get("nonpolymer_comp")
                if nonpolymer_comp:
                    chem_comp = nonpolymer_comp.get("chem_comp")
                    if chem_comp:
                        chem_id = chem_comp.get("id")
                        if chem_id:
                            chem_comp_to_pdb_map.setdefault(chem_id, set()).add(pdb_id)

    # Write to file
    output_file = "cc-to-pdb.tsv"
    with open(output_file, "w", encoding="utf-8") as f:
        for ccid, pdb_ids in chem_comp_to_pdb_map.items():
            f.write(f"{ccid}\t{' '.join(sorted(pdb_ids))}\n")

    return output_file


# Call the function
output_path = get_all_chem_comp_ids_and_write_to_file()
print(f"File saved at: {output_path}")
