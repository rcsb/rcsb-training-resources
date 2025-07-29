"""
This function downloads and cleans mmCIF data for a specific ligand in a given PDB entry.
It only retains the "chem_comp" and "atom_site" categories.

To run this script please use: python3 ligand-file-download.py
MarshalUtil is part of the RCSB PDB Python utilities package.
You can install it using: pip install rcsb.utils.io
os is a standard library for interacting with the operating system; install using pip install os
copy is used to deeply copy data containers without shared references; install using pip install copy
"""
import os
from copy import deepcopy
from rcsb.utils.io.MarshalUtil import MarshalUtil


def clean_and_collect_ligand_cif(ligand_id: str, entry_id: str, marshal_util: MarshalUtil):
    """
    This function downloads and cleans mmCIF data for a specific ligand in a given PDB entry.
    It only retains the "chem_comp" and "atom_site" categories.
    """
    # Construct URL for the mmCIF data specific to the ligand in the structure entry
    url = f"https://models.rcsb.org/v1/{entry_id}/atoms?label_comp_id={ligand_id}&encoding=cif&copy_all_categories=false&download=false"
    try:
        # Download and parse mmCIF data using the given marshal utility
        dataContainerList = marshal_util.doImport(url, fmt="mmcif")
        if not dataContainerList:
            print(f"[Warning] No data containers found for {ligand_id} in {entry_id}")
            return None

        # Extract the first data container and get all category names
        originalContainer = dataContainerList[0]
        categoryNames = originalContainer.getObjNameList()

        # Define which categories to keep
        categories_to_keep = {"chem_comp", "atom_site"}

        # Create a deep copy and remove all unwanted categories
        newContainer = deepcopy(originalContainer)
        for catName in categoryNames:
            if catName not in categories_to_keep:
                newContainer.remove(catName)

        return newContainer

    except Exception as e:
        # Print an error message if fetching or processing fails
        print(f"[Error] Failed for {ligand_id} in {entry_id}: {e}")
        return None


# This function processes a TSV input file containing ligand IDs and their associated PDB IDs.
# It collects and cleans the relevant mmCIF data for each combination, and writes it to a single output CIF file.
def process_input_file_to_single_cif(input_file_path: str, output_file: str):
    mU = MarshalUtil()
    all_containers = []

    # Open the input file and read line-by-line
    with open(input_file_path, "r", encoding="utf-8") as f:
        for line in f:
            # Each line should contain a ligand ID followed by one or more PDB entry IDs
            parts = line.strip().split()
            if len(parts) < 2:
                continue
            ligand_id = parts[0]
            pdb_ids = parts[1:]
            for entry_id in pdb_ids:
                # Fetch and clean the mmCIF data for each ligand-entry pair
                container = clean_and_collect_ligand_cif(ligand_id, entry_id, mU)
                if container:
                    # Assign a unique name to each container
                    data_name = f"{entry_id.lower()}_{ligand_id.upper()}"
                    container.setName(data_name)
                    all_containers.append(container)

    # If any cleaned containers were collected, export them to a single output CIF file
    if all_containers:
        mU.doExport(output_file, all_containers, fmt="mmcif")

        # Post-process the file to insert "#" lines after "data_" lines and remove the following line
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
                    skip_next = True  # flag the next line to be skipped
                else:
                    outfile.write(line)


# Define the input and output file paths and execute the processing function
input_file_path = "C:/Users/Krish/Documents/gitRCSB/rcsb-training-resources/cc-to-pdb.tsv"
output_file_path = "all_cleaned_ligands.cif"
process_input_file_to_single_cif(input_file_path, output_file_path)
