"""
This function queries the RCSB Search API for specific chemical component IDs (CCIDs)
and finds all PDB entries that contain each of them.
It writes a mapping from CCID to corresponding PDB IDs to a TSV file.

To run this script please use: python3 specifi-cc-to-pdb.py
These imports are from the RCSB API Python package.
You can install it using:
    pip install rcsb-api
"""
from rcsbapi.search import AttributeQuery


def get_specific_cc_to_pdb_and_write_to_file(ccids):
    """
    This function queries the RCSB Search API for specific chemical component IDs (CCIDs)
    and finds all PDB entries that contain each of them.
    It writes a mapping from CCID to corresponding PDB IDs to a TSV file.
    """
    chem_comp_to_pdb_map = {}

    # Iterate through each chemical component ID (e.g., "ATP", "GTP", etc.)
    for ccid in ccids:
        # Build a query for an exact match on the chemical component ID
        query = AttributeQuery(
            attribute="rcsb_chem_comp_container_identifiers.comp_id",  # Path to the comp_id field
            operator="exact_match",                                    # Use exact string matching
            value=ccid                                                 # Current CCID being searched
        )

        # Execute the query and collect PDB IDs that contain this CCID
        results = list(query(return_type="entry"))  # Get entries (PDB IDs) as the return type
        pdb_ids = {r for r in results}              # Convert results into a set of unique PDB IDs

        # Store the mapping: CCID → set of associated PDB IDs
        chem_comp_to_pdb_map[ccid] = pdb_ids

    # Write the final mapping to a TSV file
    # Format: <ccid>    <pdb_id1> <pdb_id2> ...
    output_file = "cc-to-pdb.tsv"
    with open(output_file, "w", encoding="utf-8") as f:
        for ccid, pdb_ids in chem_comp_to_pdb_map.items():
            f.write(f"{ccid}\t{' '.join(sorted(pdb_ids))}\n")

    return output_file


# Example usage
ccids = ["ATP", "GTP", "KYW"]  # List of CCIDs to search for
output_path = get_specific_cc_to_pdb_and_write_to_file(ccids)
print(f"File saved at: {output_path}")
