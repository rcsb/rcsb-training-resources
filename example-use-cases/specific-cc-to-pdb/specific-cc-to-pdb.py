from rcsbapi.search import AttributeQuery


def get_specific_cc_to_pdb_and_write_to_file(ccids):
    chem_comp_to_pdb_map = {}

    for ccid in ccids:
        # Query for a single CCID at a time
        query = AttributeQuery(
            attribute="rcsb_chem_comp_container_identifiers.comp_id",
            operator="exact_match",
            value=ccid
        )
        results = list(query(return_type="entry"))
        pdb_ids = {r for r in results}

        chem_comp_to_pdb_map[ccid] = pdb_ids

    # Write to TSV
    output_file = "cc-to-pdb.tsv"
    with open(output_file, "w", encoding="utf-8") as f:
        for ccid, pdb_ids in chem_comp_to_pdb_map.items():
            f.write(f"{ccid}\t{' '.join(sorted(pdb_ids))}\n")

    return output_file


# Example usage
ccids = ["ATP", "GTP", "KYW"]
output_path = get_specific_cc_to_pdb_and_write_to_file(ccids)
print(f"File saved at: {output_path}")
