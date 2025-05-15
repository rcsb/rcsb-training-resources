"""
Extract protein structures with small molecules bound at transmembrane regions from the RCSB PDB.

General workflow:
    1. First, use the search API to get a list of all entity instances (chains) of transmembrane proteins (based on PDBTM annotations)
    2. Then, use the data API to fetch all feature data for those instances (chains) (“rcsb_polymer_instance_feature”)
    3. Next, filter down the feature data to only retain the following:
        - All PDBTM features (which gives you the sequence positions of TM segments)
        - All ligand interactions (along with the sequence positions of those interactions)
    4. Last, write custom code to parse these data to identify which ligands interactions overlap with TM segments

To run:
    python3 transmembrane_protein_ligand_interactions.py

"""

from rcsbapi.search import search_attributes as attrs
from rcsbapi.data import DataQuery as Query
from pprint import pprint


def main():
    ## 1. First, use the search API to get a list of all entity instances (chains) with PDBTM annotations (https://pdbtm.unitmp.org/)
    q1 = attrs.rcsb_polymer_entity_annotation.type == "PDBTM"
    search_results = list(q1(return_type="polymer_instance"))
    print(f"# Entity instances: {len(search_results)}")
    print(f"First few results: {search_results[:10]}")  # Print the first 10 resulting IDs


    ## 2. Then, use the data API to fetch all feature data for those instances (chains) (“rcsb_polymer_instance_feature”)
    ##    This might take a while if the length of input_ids is very large; consider testing with just the first 10 IDs.
    data_query = Query(
        input_type="polymer_entity_instances",
        input_ids=search_results,  # For testing, only use a few input IDs, e.g., search_results[:10]
        return_data_list=["rcsb_polymer_instance_feature", "polymer_entity_instances.rcsb_id"]
    )
    data_results = data_query.exec()
    # pprint(data_results["data"]["polymer_entity_instances"][0])  # Print the first result


    ## 3. Next, filter down the feature data to only retain the the PDBTM and ligand interaction information
    instance_feature_data_filtered = {}  # Will be a dictionary with instance IDs as keys and filtered feature data as values
    for instance in data_results["data"]["polymer_entity_instances"]:
        instance_id = instance["rcsb_id"]
        filtered_feature_data = {"ligand_interactions": [], "membrane_segments": []}
        try:
            feature_list = instance.get("rcsb_polymer_instance_feature")
            if feature_list:
                for feature in feature_list:
                    if feature["type"] == "MEMBRANE_SEGMENT" and feature["provenance_source"] == "PDBTM":
                        filtered_feature_data["membrane_segments"].append(extract_feature_data(feature))
                    if feature["type"] == "LIGAND_INTERACTION":
                        filtered_feature_data["ligand_interactions"].append(extract_feature_data(feature))
                instance_feature_data_filtered.update({instance_id: filtered_feature_data})
        except Exception as e:
            print(f"Failing for {instance_id} with: {e}")
    #
    pprint(list(instance_feature_data_filtered.items())[0])  # print out the first filtered result


    ## 4. Last, write custom code to parse the filtered feature data generated above for your specific research task,
    ##    such as to identify which ligand interactions overlap with TM segments.
    ##    The data returned by the above step will have the following structure (using 1BY3.A as an example):

    # >>> print(instance_feature_data_filtered)
    # {
    #     '1BY3.A': {
    #         'ligand_interactions': [
    #             {
    #                 'description': 'Software generated binding site for ligand entity 2 component OES instance B',
    #                 'provenance_source': 'PDB',
    #                 'additional_properties': [{'values': ['B'], 'name': 'PARTNER_ASYM_ID'}, {'values': ['OES'], 'name': 'PARTNER_COMP_ID'}],
    #                 'name': 'ligand OES',
    #                 'feature_positions': [
    #                     {'beg_seq_id': 200, 'beg_comp_id': 'SER'},
    #                     {'beg_seq_id': 176, 'beg_comp_id': 'THR'},
    #                     {'beg_seq_id': 211, 'beg_comp_id': 'GLN'},
    #                     {'beg_seq_id': 174, 'beg_comp_id': 'PHE'},
    #                     {'beg_seq_id': 198, 'beg_comp_id': 'ALA'},
    #                     {'beg_seq_id': 213, 'beg_comp_id': 'TYR'},
    #                     {'beg_seq_id': 209, 'beg_comp_id': 'GLU'}
    #                 ]
    #             },
    #             {
    #                 'description': 'Software generated binding site for ligand entity 2 component OES instance C',
    #                 'provenance_source': 'PDB',
    #                 'additional_properties': [{'values': ['C'], 'name': 'PARTNER_ASYM_ID'}, {'values': ['OES'], 'name': 'PARTNER_COMP_ID'}],
    #                 'name': 'ligand OES',
    #                 'feature_positions': [
    #                     {'beg_seq_id': 676, 'beg_comp_id': 'VAL'},
    #                     {'beg_seq_id': 710, 'beg_comp_id': 'ALA'},
    #                     {'beg_seq_id': 168, 'beg_comp_id': 'ALA'},
    #                     {'beg_seq_id': 166, 'beg_comp_id': 'PHE'},
    #                     {'beg_seq_id': 708, 'beg_comp_id': 'ALA'},
    #                     {'beg_seq_id': 174, 'beg_comp_id': 'PHE'},
    #                     {'beg_seq_id': 175, 'beg_comp_id': 'GLN'},
    #                     {'beg_seq_id': 167, 'beg_comp_id': 'LYS'}
    #                 ]
    #             },
    #             ...
    #         ],
    #         'membrane_segments': [
    #             {
    #                 'provenance_source': 'PDBTM',
    #                 'feature_positions': [
    #                     {'beg_seq_id': 28, 'end_seq_id': 75},
    #                     {'beg_seq_id': 120, 'end_seq_id': 140},
    #                     {'beg_seq_id': 146, 'end_seq_id': 155},
    #                     {'beg_seq_id': 161, 'end_seq_id': 168},
    #                     {'beg_seq_id': 172, 'end_seq_id': 180},
    #                     {'beg_seq_id': 191, 'end_seq_id': 198},
    #                     {'beg_seq_id': 212, 'end_seq_id': 220},
    #                     {'beg_seq_id': 227, 'end_seq_id': 235},
    #                     {'beg_seq_id': 279, 'end_seq_id': 286},
    #                     {'beg_seq_id': 294, 'end_seq_id': 301},
    #                     {'beg_seq_id': 355, 'end_seq_id': 362},
    #                     {'beg_seq_id': 372, 'end_seq_id': 379},
    #                     {'beg_seq_id': 432, 'end_seq_id': 439},
    #                     {'beg_seq_id': 445, 'end_seq_id': 452},
    #                     {'beg_seq_id': 476, 'end_seq_id': 483},
    #                     {'beg_seq_id': 489, 'end_seq_id': 496},
    #                     {'beg_seq_id': 519, 'end_seq_id': 526},
    #                     {'beg_seq_id': 534, 'end_seq_id': 540},
    #                     {'beg_seq_id': 568, 'end_seq_id': 575},
    #                     {'beg_seq_id': 582, 'end_seq_id': 589},
    #                     {'beg_seq_id': 612, 'end_seq_id': 619},
    #                     {'beg_seq_id': 630, 'end_seq_id': 637},
    #                     {'beg_seq_id': 656, 'end_seq_id': 663},
    #                     {'beg_seq_id': 674, 'end_seq_id': 682},
    #                     {'beg_seq_id': 704, 'end_seq_id': 712}
    #                 ]
    #             }
    #         ]
    #     }
    # }


def extract_feature_data(featureD):
    """Extract out relevant information from ligand interaction or membrane segment features"""
    if featureD["type"] == "LIGAND_INTERACTION":
        feature_keys = ["description", "provenance_source", "additional_properties", "name", "feature_positions"]
        feature_position_keys = ["beg_seq_id", "beg_comp_id"]
    elif featureD["type"] == "MEMBRANE_SEGMENT":
        feature_keys = ["provenance_source", "feature_positions"]
        feature_position_keys = ["beg_seq_id", "end_seq_id"]
    else:
        raise ValueError("Set of feature information to extract not defined for feature 'type' %s" % featureD["type"])
    #
    extracted_feature_data = {}
    for k in feature_keys:
        feature_data = featureD.get(k, None)
        if k == "feature_positions":
            feature_positions = []
            for position in feature_data:
                feature_positions.append({fpk: fpv for fpk, fpv in position.items() if fpk in feature_position_keys})
            extracted_feature = feature_positions
        else:
            extracted_feature = feature_data
        #
        extracted_feature_data[k] = extracted_feature
    #
    return extracted_feature_data


if __name__ == "__main__":
    main()
