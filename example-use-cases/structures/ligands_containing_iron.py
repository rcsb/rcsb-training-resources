"""
Search for proteins containing iron ions as ligands (as free ion or bound to a heme or iron-sulfur cluster etc.)
and return non-redundant structures with the highest resolution.

General workflow:
    1. First, use the search API to get a list of all chemical components containing iron
    2. Then, use the search API to find PDB sequences that have components of interest as ligands in proximity
    3. Next, remove sequences that share at least 95% identity and return sequences from structures with the
       highest resolution

To run:
    python3 ligands_containing_iron.py

"""

import time
from rcsbapi.search import AttributeQuery, ChemSimilarityQuery, GroupBy, RankingCriteriaType


def search_iron_containing_ccd():
    chem_query = ChemSimilarityQuery(
        value="Fe",
        query_type="formula",
        match_subset=True
    )
    return list(chem_query(return_type="mol_definition", results_verbosity="compact"))


def search_sequence_targets(ligands):

    q1 = AttributeQuery(
        attribute="rcsb_ligand_neighbors.ligand_comp_id",
        operator="in",
        value=ligands
    )

    q2 = AttributeQuery(
        attribute="entity_poly.rcsb_entity_polymer_type",
        operator="exact_match",
        value="Protein"
    )

    q3 = AttributeQuery(
        attribute="exptl.method",
        operator="exact_match",
        value="X-RAY DIFFRACTION"
    )

    group_by = GroupBy(
        aggregation_method="sequence_identity",
        similarity_cutoff=95,
        ranking_criteria_type=RankingCriteriaType(
            sort_by="rcsb_entry_info.resolution_combined",
            direction="asc"
        )
    )

    query = q1 & q2 & q3
    return list(query(
        return_type="polymer_entity",
        group_by=group_by,
        group_by_return_type="representatives"
    ))


if __name__ == "__main__":
    s_time = time.time()
    ccd_ids = search_iron_containing_ccd()
    final_sequences = search_sequence_targets(ccd_ids)
    print(f"Total number of results: {len(final_sequences)}")
    print(f"List of results: {final_sequences}")
    e_time = time.time()
    print(f"Total execution time: {e_time - s_time:.4f} seconds")
