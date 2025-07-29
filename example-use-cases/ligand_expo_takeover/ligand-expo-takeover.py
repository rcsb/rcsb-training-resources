"""
This function queries all PDB structures to extract chemical component IDs (ligands)
from branched, polymer, and nonpolymer entities, and writes a mapping of
chem_comp_id -> associated PDB IDs to a TSV file.

To run this script please use: python3 ligand-expo-takeover.py
These imports are from the RCSB API Python package.
You can install it using: pip install rcsb-api
"""
from rcsbapi.data import DataQuery as Query
from rcsbapi.data import ALL_STRUCTURES


def get_all_chem_comp_ids_and_write_to_file():
    """
    This function queries all PDB structures to extract chemical component IDs (ligands)
    from branched, polymer, and nonpolymer entities, and writes a mapping of
    chem_comp_id -> associated PDB IDs to a TSV file.
    """
    # Initialize the data query to retrieve relevant chemical component data
    query = Query(
        input_type="entries",              # Query all structure entries
        input_ids=ALL_STRUCTURES,         # Constant representing all known structures
        return_data_list=[
            "rcsb_id",  # PDB ID of the structure
            "polymer_entities.chem_comp_nstd_monomers.chem_comp.id",   # Non-standard monomers in polymer chains
            "branched_entities.chem_comp_monomers.chem_comp.id",       # Monomers in branched entities
            "nonpolymer_entities.nonpolymer_comp.chem_comp.id"         # Ligands in nonpolymer entities
        ]
    )

    # Execute the query with a progress bar
    result = query.exec(progress_bar=True)

    # Extract list of returned structure entries
    entries = result.get("data", {}).get("entries", [])

    # Dictionary to collect mapping from chem_comp_id to a set of PDB IDs
    chem_comp_to_pdb_map = {}

    # Iterate over all returned entries
    for entry in entries:
        pdb_id = entry.get("rcsb_id")

        # --- 1. Branched Entities ---
        # These may include things like saccharides
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
                                # Add chem_comp_id → pdb_id to the mapping
                                chem_comp_to_pdb_map.setdefault(chem_id, set()).add(pdb_id)

        # --- 2. Polymer Entities ---
        # Includes non-standard residues within polymer chains
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

        # --- 3. Nonpolymer Entities ---
        # Usually small molecule ligands
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

    # Write the final mapping to a TSV file
    # Format: <chem_comp_id>    <pdb_id1> <pdb_id2> ...
    output_file = "cc-to-pdb.tsv"
    with open(output_file, "w", encoding="utf-8") as f:
        for ccid, pdb_ids in chem_comp_to_pdb_map.items():
            f.write(f"{ccid}\t{' '.join(sorted(pdb_ids))}\n")

    return output_file


# Run the function and report the output file path
output_path = get_all_chem_comp_ids_and_write_to_file()
print(f"File saved at: {output_path}")
