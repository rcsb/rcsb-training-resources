"""
This script extracts the coordinate data for a given list of chemical component (ligand) IDs
as present in PDB entry files. The extracted mmCIF data is cleaned to retain only the
"chem_comp" and "atom_site" categories.

Requirements:
    pip install "rcsb-api>=1.4.0"
    pip install mmcif

Usage:
    # Get usage details
        python3 extract_ligand_coordinates.py --help

    # Extract the HEM ligand coordinates from a limit of 10 PDB entries
        python3 extract_ligand_coordinates.py -c HEM -n 10

    # Extract the HEM and PO4 ligand coordinates from PDB entries 4HHB and 3HHB,
    # and save output to "output_files" directory
        python3 extract_ligand_coordinates.py -c HEM PO4 -p 4HHB 3HHB -o output_files

    # Extract the CPT ligand coordinates from all PDB entries that contain it
        python3 extract_ligand_coordinates.py -c CPT

Output:
    output/<CID>-coordinates.cif  # one file per CID
"""

from pathlib import Path
import argparse
import io
import time
import logging
from mmcif.io.PdbxReader import PdbxReader
from mmcif.io.PdbxWriter import PdbxWriter
from rcsbapi.search import AttributeQuery
from rcsbapi.model import ModelQuery


logging.basicConfig(level=logging.WARNING, format="%(asctime)s [%(levelname)s]-%(module)s.%(funcName)s: %(message)s")
logger = logging.getLogger(__name__)


class LigandCoordinatesExtract:
    """
    This class provides methods to download and clean mmCIF data for specific ligands in PDB entries.
    It retains only the "chem_comp" and "atom_site" categories.
    """

    def __init__(self):
        self.model_query = ModelQuery()
        self.data_container_list = []

    def process_one(self, ccd_id: str, pdb_id: str):
        """
        Retrieve and cleans mmCIF data for a specific ligand in a given PDB entry.
        """
        try:
            # Query the API for atoms related to the ligand in the specified entry
            result = self.model_query.get_atoms(entry_id=pdb_id, label_comp_id=ccd_id)
            if not result:
                logger.warning(f"No data found for {ccd_id} in {pdb_id}")
                return None

            # Read the result into a file-like object
            file_like = io.StringIO(result)
            pR = PdbxReader(file_like)
            pR.read(self.data_container_list, ["chem_comp", "atom_site"])
        except Exception as e:
            logger.error(f"Failed for {ccd_id} in {pdb_id}: {e}")
            return None

    def process_all(self, ccd_id: str, ccd_pdbid_list: list, fp_out: str, write_interval: int = 10):
        """
        Processes a dictionary mapping ligand IDs to lists of PDB IDs.
        Collects and cleans the relevant mmCIF data for each combination.
        """
        write_i = 0
        append_flag = False
        for pdb_id in ccd_pdbid_list:
            logger.debug(f"Retrieving {ccd_id} coordinates from {pdb_id}")
            self.process_one(ccd_id, pdb_id)
            write_i += 1
            if write_i == write_interval:
                self.write_data(fp_out, append_flag)
                write_i = 0
                append_flag = True
        if len(self.data_container_list) > 0:
            self.write_data(fp_out, append_flag)

    def write_data(self, fp_out: str, append_flag: bool = True):
        """
        Writes the cleaned mmCIF data to the specified output file.
        """
        try:
            write_mode = "a" if append_flag else "w"
            if append_flag:
                write_mode = "a"
                print(f"Appending output to '{fp_out}' for PDB IDs {[dc.getName() for dc in self.data_container_list]}")
            else:
                write_mode = "w"
                print(f"Writing output to '{fp_out}' for PDB IDs {[dc.getName() for dc in self.data_container_list]}")
            #
            with open(fp_out, write_mode) as f:
                pW = PdbxWriter(f)
                pW.write(self.data_container_list)
            #
            self.data_container_list = []
        #
        except Exception as e:
            logger.error(f"Failed to write output with exception: {e}")


def search_one_cid(ccd_id):
    """Retrieve all PDB IDs with a specific CCD ID

    Args:
        ccd_id (str): CCD ID (e.g., ATP)

    Returns:
        list: PDB IDs list
    """
    # Build a query for an exact match on the chemical component ID
    query = AttributeQuery(
        attribute="rcsb_chem_comp_container_identifiers.comp_id",  # Path to the comp_id field
        operator="exact_match",                                    # Use exact string matching
        value=ccd_id,
    )

    # Execute the query and collect PDB IDs that contain this CCID
    results = list(query(return_type="entry"))  # Get entries (PDB IDs) as the return type
    return [pdb_id.upper() for pdb_id in results]  # Convert results into a set of unique PDB IDs


def search_pdb_by_ccid(ccd_id, pdbid_limit_list=None, pdb_limit_num=None):
    """search PDB IDs with the CCD ID

    Args:
        ccd_id (str): CCD ID to search
        pdbid_limit_list (list): list of PDB IDs to limit extraction to
        pdb_limit_num (int): max number of PDB IDs to extract the ligand from

    Returns:
        l_ccd_pdbids: list of PDB IDs
    """
    l_ccd_pdbids = []
    l_pdb = search_one_cid(ccd_id)
    if pdbid_limit_list:
        for pdb_id in pdbid_limit_list:
            if pdb_id.upper() in l_pdb:
                l_ccd_pdbids.append(pdb_id.upper())
        print(f"Will extract {ccd_id} coordinates from PDB IDs based on input IDs ({len(l_ccd_pdbids)} total)")
    elif pdb_limit_num:
        if len(l_pdb) <= pdb_limit_num:
            l_ccd_pdbids = l_pdb
        else:
            l_ccd_pdbids = l_pdb[:pdb_limit_num]
        print(f"Will extract {ccd_id} coordinates from PDB IDs based on input limit ({len(l_ccd_pdbids)} total)")
    else:
        l_ccd_pdbids = l_pdb
        print(f"Will extract {ccd_id} coordinates from all CCID-containing PDB IDs ({len(l_ccd_pdbids)} total)")
    #
    if len(l_ccd_pdbids) >= 25:
        print(f"  PDB IDs (first 25): {l_ccd_pdbids[:25]}")
    else:
        print(f"  PDB IDs: {l_ccd_pdbids[:25]}")

    return l_ccd_pdbids


def extract_ligand_coordinates(ccid_list: list, pdbid_limit_list: list, pdb_limit_num: int, output_dir: str, write_interval: int = 10):
    """Extract ligand coordinates from PDB entries for provided list of CCD IDs and specified limits.

    Args:
        ccid_list (list): list of CC IDs to extract
        pdbid_limit_list (list): list of PDB IDs to limit extraction to
        pdb_limit_num (int): max number of PDB IDs to extract the ligand from
        output_dir (str): output directory to write extracted ligand coordinate files to
        write_interval (int): how many sets of coordinates to write out at once (and clear internal memory)
    """
    # Create the output directory if it doesn't exist
    Path(output_dir).mkdir(parents=True, exist_ok=True)
    print(f"Will be writing output files to directory: {output_dir}")

    for ccd_id in ccid_list:
        fp_out = Path(output_dir) / f"{ccd_id}-coordinates.cif"
        ccd_pdbid_list = search_pdb_by_ccid(ccd_id, pdbid_limit_list, pdb_limit_num)
        if ccd_pdbid_list:
            ligand_extractor = LigandCoordinatesExtract()
            ligand_extractor.process_all(ccd_id, ccd_pdbid_list, fp_out, write_interval=write_interval)
        else:
            logger.error(f"Failed to find PDB entries with CCD ID {ccd_id}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Extract coordinates of specific ligands from PDB entries")

    # take mandatory author-provided CCD IDs, can process multiple IDs
    parser.add_argument("-c", "--ccids", nargs="+", required=True, help="Space-separated list of ligand CCD IDs (e.g. CPT)")
    parser.add_argument("-o", "--output-dir", default="output", required=False, help="Path/directory to save output files")

    # add mutually exclusive options limitating extraction by PDB IDs or number of entries
    limit_extraction = parser.add_mutually_exclusive_group(required=False)
    limit_extraction.add_argument("-p", "--pdbids", nargs="+", default=None, help="Space-separated list of PDB IDs to limit the extraction on")
    limit_extraction.add_argument("-n", "--max-num-ids", default=None, help="Limit extraction to the specified number of PDB entries")

    args = parser.parse_args()

    input_ccid_list = [ccid.upper() for ccid in args.ccids]
    input_pdbid_list = [pdbid.upper() for pdbid in args.pdbids] if args.pdbids else None
    n_pdb_limit = int(args.max_num_ids) if args.max_num_ids else None
    output_dir = args.output_dir
    print(f"List of CCD IDs for which to extract ligand coordinates: {','.join(input_ccid_list)}")

    start = time.time()
    extract_ligand_coordinates(input_ccid_list, input_pdbid_list, n_pdb_limit, output_dir)
    end = time.time()
    print(f"Processing completed in {end - start:.2f} seconds.")
