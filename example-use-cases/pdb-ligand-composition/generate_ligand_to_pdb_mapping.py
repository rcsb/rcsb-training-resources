"""
This script generates a mapping between all ligands in the PDB and the corresponding
list of PDB structures in which they are present. The complete mapping is written to
an output TSV file.

By default, the list of corresponding PDB IDs includes those in which the ligand
(or, "chemical component") exists in one of the following forms:

    - As a nonpolymer entity (e.g., a non-covalently bound ATP ligand)
    - As a monomer within a branched polymer entity (e.g., a subunit in a carbohydrate)
    - As a non-standard amino acid monomer within a polymer entity

To exclude any of these from the mapping, comment out the relevant items in the
assignment of `chemical_component_types_to_include` below.


Requirements:
    pip install "rcsb-api>=1.4.0"

Usage:
    python3 generate_ligand_to_pdb_mapping.py

Output:
    cc-to-pdb.tsv  # Or custom name as can be specified below

    # Format: <chem_comp_id>    <pdb_id1> <pdb_id2> ...
"""

import time
from rcsbapi.config import config
from rcsbapi.data import DataQuery as Query
from rcsbapi.data import ALL_STRUCTURES


# Output and runtime settings - configure as desired
output_file = "cc-to-pdb.tsv"

chemical_component_types_to_include = [
    "nonpolymer_entities.rcsb_nonpolymer_entity_container_identifiers.nonpolymer_comp_id",  # Ligands in nonpolymer entities
    "branched_entities.rcsb_branched_entity_container_identifiers.chem_comp_monomers",      # Monomers in branched entities
    "polymer_entities.rcsb_polymer_entity_container_identifiers.chem_comp_nstd_monomers",   # Non-standard monomers in polymer chains
]

config.DATA_API_MAX_CONCURRENT_REQUESTS = 10  # Specify custom API settings to speed up fetch


def fetch_all_chem_comp_ids(chem_comp_types_to_include: list):
    """Fetch the chemical component and PDB mapping data from RCSB.org using the Data API"""

    # Initialize the data query to retrieve relevant chemical component data
    query = Query(
        input_type="entries",             # Query all structure entries
        input_ids=ALL_STRUCTURES,         # Constant representing all known structures
        return_data_list=["rcsb_id"]+chem_comp_types_to_include
    )

    # Execute the query with a progress bar
    result = query.exec(progress_bar=True)

    # Extract list of returned structure entries
    entry_chem_comp_results = result.get("data", {}).get("entries", [])

    return entry_chem_comp_results


def process_chem_comp_results_and_write_to_file(entry_chem_comp_results: dict, output_file: str):
    """Process the fetched Data API data to create the mapping between chemical component IDs and PDB IDs"""

    # Dictionary to collect mapping from chem_comp_id to a set of PDB IDs
    chem_comp_to_pdb_map = {}

    # Iterate over all returned entries
    for entry in entry_chem_comp_results:
        pdb_id = entry.get("rcsb_id")

        # --- 1. Nonpolymer Entities ---
        # Small molecule ligands
        nonpolymer_entities = entry.get("nonpolymer_entities")
        if nonpolymer_entities:
            for nonpolymer in nonpolymer_entities:
                nonpolymer_comp = nonpolymer.get("rcsb_nonpolymer_entity_container_identifiers")
                if nonpolymer_comp:
                    chem_id = nonpolymer_comp.get("nonpolymer_comp_id")
                    if chem_id:
                        chem_comp_to_pdb_map.setdefault(chem_id, set()).add(pdb_id)

        # --- 2. Branched Entities ---
        # Includes things like monomers in saccharides
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

        # --- 3. Polymer Entities ---
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

    # Write the final mapping to a TSV file
    with open(output_file, "w", encoding="utf-8") as f:
        for ccid, pdb_ids in chem_comp_to_pdb_map.items():
            f.write(f"{ccid}\t{' '.join(sorted(pdb_ids))}\n")

    print(f"File saved at: {output_file}")


if __name__ == "__main__":
    start = time.time()
    entry_chem_comp_data = fetch_all_chem_comp_ids(chemical_component_types_to_include)
    process_chem_comp_results_and_write_to_file(entry_chem_comp_data, output_file)
    end = time.time()
    print(f"Processing completed in {end - start:.2f} seconds.")
