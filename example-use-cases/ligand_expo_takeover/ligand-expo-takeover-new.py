"""
This function queries all PDB structures to extract chemical component IDs (ligands)
from branched, polymer, and nonpolymer entities, and writes a mapping of
chem_comp_id -> associated PDB IDs to a TSV file.

To run this script please use: python3 ligand-expo-takeover.py
These imports are from the RCSB API Python package.
You can install it using: pip install rcsb-api


*** CHENGHUA: You should install it as follows (for now): ***
pip install "git+https://github.com/rcsb/py-rcsb-api.git@dev-httpx#egg=rcsb-api"

"""
from rcsbapi.config import config
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
            "polymer_entities.rcsb_polymer_entity_container_identifiers.chem_comp_nstd_monomers",   # Non-standard monomers in polymer chains
            "branched_entities.rcsb_branched_entity_container_identifiers.chem_comp_monomers",       # Monomers in branched entities
            "nonpolymer_entities.rcsb_nonpolymer_entity_container_identifiers.nonpolymer_comp_id"         # Ligands in nonpolymer entities
            # Below are slower to fetch than above
            # "polymer_entities.chem_comp_nstd_monomers.chem_comp.id",   # Non-standard monomers in polymer chains
            # "branched_entities.chem_comp_monomers.chem_comp.id",       # Monomers in branched entities
            # "nonpolymer_entities.nonpolymer_comp.chem_comp.id"         # Ligands in nonpolymer entities
        ]
    )

    # ** CHENGHUA **
    # Specify custom configuration settings to speed up fetch
    config.DATA_API_BATCH_ID_SIZE = 200
    config.DATA_API_MAX_CONCURRENT_REQUESTS = 10

    # Execute the query with a progress bar
    result = query.exec(progress_bar=True)

    # Extract list of returned structure entries
    entries = result.get("data", {}).get("entries", [])

    # Dictionary to collect mapping from chem_comp_id to a set of PDB IDs
    chem_comp_to_pdb_map = {}

    # Iterate over all returned entries
    for entry in entries:
        pdb_id = entry.get("rcsb_id")

        # --- 1. Polymer Entities ---
        # Includes non-standard residues within polymer chains
        polymer_entities = entry.get("polymer_entities")
        if polymer_entities:
            for polymer in polymer_entities:
                polymer_entity_container_identifiers = polymer.get("rcsb_polymer_entity_container_identifiers")
                if polymer_entity_container_identifiers:
                    chem_comp_nstd_monomers = polymer_entity_container_identifiers.get("chem_comp_nstd_monomers")
                    if chem_comp_nstd_monomers:
                        for chem_id in chem_comp_nstd_monomers:
                            # Add chem_comp_id → pdb_id to the mapping
                            chem_comp_to_pdb_map.setdefault(chem_id, set()).add(pdb_id)

        # --- 2. Branched Entities ---
        # These may include things like saccharides
        branched_entities = entry.get("branched_entities")
        if branched_entities:
            for branched in branched_entities:
                branched_container_identifiers = branched.get("rcsb_branched_entity_container_identifiers")
                if branched_container_identifiers:
                    chem_comp_monomers = branched_container_identifiers.get("chem_comp_monomers")
                    if chem_comp_monomers:
                        for chem_id in chem_comp_monomers:
                            # Add chem_comp_id → pdb_id to the mapping
                            chem_comp_to_pdb_map.setdefault(chem_id, set()).add(pdb_id)

        # --- 3. Nonpolymer Entities ---
        # Usually small molecule ligands
        nonpolymer_entities = entry.get("nonpolymer_entities")
        if nonpolymer_entities:
            for nonpolymer in nonpolymer_entities:
                nonpolymer_comp = nonpolymer.get("rcsb_nonpolymer_entity_container_identifiers")
                if nonpolymer_comp:
                    chem_id = nonpolymer_comp.get("nonpolymer_comp_id")
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
