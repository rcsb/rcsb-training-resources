"""
pip install rcsb-api
"""
import csv
from rcsbapi.search import StructMotifQuery, StructMotifResidue
from rcsbapi.data import DataQuery

EC_CLASS_MAP = {
    "1": "Oxidoreductase",
    "2": "Transferase",
    "3": "Hydrolase",
    "4": "Lyase",
    "5": "Isomerase",
    "6": "Ligase",
    "7": "Translocase",
}

# Pick EC numbers, preferring UniProt
def pick_ec(ec_data, ec_list):
    if not ec_data:
        return
    
    uniprot_ecs = []
    fallback_ecs = []

    for item in ec_data:
        ec = item.get("ec")
        if not ec:
            continue
        if item.get("provenance_source") == "UniProt":
            uniprot_ecs.append(ec)
        else:
            fallback_ecs.append(ec)

    if uniprot_ecs:
        ec_list.update(uniprot_ecs)
    else:
        ec_list.update(fallback_ecs)

# Convert EC numbers to unique enzyme class names
def ec_to_classes(ec_list):
    classes = set()

    for ec in ec_list:
        if not ec:
            continue
        first_digit = ec.split(".", 1)[0]
        enzyme_class = EC_CLASS_MAP.get(first_digit)
        if enzyme_class:
            classes.add(enzyme_class)

    return sorted(classes)

if __name__ == "__main__":
    res1 = StructMotifResidue(
        struct_oper_id="1",
        chain_id="A",
        label_seq_id=41
    )
    res2 = StructMotifResidue(
        struct_oper_id="1",
        chain_id="A",
        label_seq_id=84
    )
    res3 = StructMotifResidue(
        struct_oper_id="1",
        chain_id="A",
        label_seq_id=177
    )
    res4 = StructMotifResidue(
        struct_oper_id="1",
        chain_id="A",
        label_seq_id=178
    )
    res5 = StructMotifResidue(
        struct_oper_id="1",
        chain_id="A",
        label_seq_id=179
    )
    res6 = StructMotifResidue(
        struct_oper_id="1",
        chain_id="A",
        label_seq_id=180
    )
    residues = [res1, res2, res3, res4, res5, res6]
    search_query = StructMotifQuery(
        entry_id="1PQ5",
        residue_ids=residues,
        atom_pairing_scheme="ALL",
        rmsd_cutoff=2
    )
    results=list(search_query(
        results_verbosity="verbose",
        return_type="assembly"
    ))
    print(f"Total number of assemblies: {len(results)}")

    all_rows = []

    for id in results:
        entry_id = id['identifier'].split("-")[0].lower()
        extended_id = f"pdb_{entry_id.zfill(8)}"
        context = id['services'][0]['nodes'][0]['match_context']

        for c in context:
            score = c['score']

            # unique chain IDs → string
            residue_ids = c['residue_ids']
            label_asym_ids = sorted({r['label_asym_id'] for r in residue_ids})
            chain_ids = ",".join(label_asym_ids)

            # Prepare sets
            protein = set()
            ec = set()

            # Fetch additional metadata per chain
            for asym_id in label_asym_ids:
                data_query = DataQuery(
                    input_type="polymer_entity_instance",
                    input_ids={"entry_id": entry_id, "asym_id": asym_id},
                    return_data_list=[
                        "polymer_entity.rcsb_polymer_entity_container_identifiers.entity_id",
                        "polymer_entity.rcsb_polymer_entity_container_identifiers.uniprot_ids",
                        "polymer_entity.rcsb_polymer_entity.pdbx_description",
                        "polymer_entity.rcsb_polymer_entity.rcsb_enzyme_class_combined.ec",
                        "polymer_entity.rcsb_polymer_entity.rcsb_enzyme_class_combined.provenance_source"
                    ]
                )
                data = data_query.exec()
                entity = data['data']['polymer_entity_instances'][0]['polymer_entity']

                entity_id = entity['rcsb_polymer_entity_container_identifiers']['entity_id']
                uniprot_ids = entity['rcsb_polymer_entity_container_identifiers']['uniprot_ids']
                rcsb_polymer_entity = entity['rcsb_polymer_entity']

                if rcsb_polymer_entity.get('pdbx_description'):
                    protein.add(rcsb_polymer_entity['pdbx_description'])

                if rcsb_polymer_entity.get('rcsb_enzyme_class_combined'):
                    pick_ec(rcsb_polymer_entity['rcsb_enzyme_class_combined'], ec)

            # Append row with helper key for deduplication
            all_rows.append({
                "_key": (entry_id, entity_id),
                "PDB ID": extended_id,
                "Entyty ID": entity_id,
                "Chain ID(s)": chain_ids,
                "UniProt ID(s)": ", ".join(sorted(uniprot_ids)) if uniprot_ids else "",
                "Protein(s)": ", ".join(sorted(protein)) if protein else "",
                "EC Class": ", ".join(ec_to_classes(ec)) if ec else "",
                "EC Number(s)": ", ".join(sorted(ec)) if ec else "",
                "RMSD": score
            })

    # Keep only lowest RMSD per (entry_id, entity_id)
    best_rows = {}
    for row in all_rows:
        key = row["_key"]
        if key not in best_rows or row["RMSD"] < best_rows[key]["RMSD"]:
            best_rows[key] = row

    # Prepare final rows
    final_rows = []
    for row in best_rows.values():
        row.pop("_key", None)  # remove helper key
        final_rows.append(row)

    # Write CSV
    with open("/Users/yana.rose/trypsin-catalytic-site.csv", "w") as handle:
        headers = list(final_rows[0].keys())
        writer = csv.DictWriter(handle, fieldnames=headers)
        writer.writeheader()
        writer.writerows(final_rows)

    print("Wrote data to trypsin-catalytic-site.csv")