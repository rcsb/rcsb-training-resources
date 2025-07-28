import os
from copy import deepcopy
from rcsb.utils.io.MarshalUtil import MarshalUtil


def clean_and_collect_ligand_cif(ligand_id: str, entry_id: str, marshal_util: MarshalUtil):
    url = f"https://models.rcsb.org/v1/{entry_id}/atoms?label_comp_id={ligand_id}&encoding=cif&copy_all_categories=false&download=false"
    try:
        dataContainerList = marshal_util.doImport(url, fmt="mmcif")
        if not dataContainerList:
            print(f"[Warning] No data containers found for {ligand_id} in {entry_id}")
            return None

        originalContainer = dataContainerList[0]
        categoryNames = originalContainer.getObjNameList()
        categories_to_keep = {"chem_comp", "atom_site"}

        newContainer = deepcopy(originalContainer)
        for catName in categoryNames:
            if catName not in categories_to_keep:
                newContainer.remove(catName)

        return newContainer

    except Exception as e:
        print(f"[Error] Failed for {ligand_id} in {entry_id}: {e}")
        return None


def process_input_file_to_single_cif(input_file_path: str, output_file: str):
    mU = MarshalUtil()
    all_containers = []

    with open(input_file_path, "r", encoding="utf-8") as f:
        for line in f:
            parts = line.strip().split()
            if len(parts) < 2:
                continue
            ligand_id = parts[0]
            pdb_ids = parts[1:]
            for entry_id in pdb_ids:
                container = clean_and_collect_ligand_cif(ligand_id, entry_id, mU)
                if container:
                    data_name = f"{entry_id.lower()}_{ligand_id.upper()}"
                    container.setName(data_name)
                    all_containers.append(container)

    if all_containers:
        mU.doExport(output_file, all_containers, fmt="mmcif")
        with open(output_file, "r", encoding="utf-8") as infile:
            lines = infile.readlines()

        with open(output_file, "w", encoding="utf-8") as outfile:
            skip_next = False
            for i, line in enumerate(lines):
                if skip_next:
                    skip_next = False
                    continue  # skip the line that came after "data_"

                if line.strip().startswith("data_"):
                    outfile.write(line)
                    outfile.write("#\n")
                    skip_next = True  # remove the line immediately after
                else:
                    outfile.write(line)


input_file_path = "ex.tsv"
output_file_path = "all_cleaned_ligands.cif"
process_input_file_to_single_cif(input_file_path, output_file_path)
