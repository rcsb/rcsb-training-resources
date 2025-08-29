"""
This function downloads and cleans mmCIF data for a specific ligand in a given PDB entry.
It only retains the "chem_comp" and "atom_site" categories.

To run this script please use: python3 ligand-file-download-api.py
"""
import argparse
import io
import time
from mmcif.io.PdbxReader import PdbxReader
from mmcif.io.PdbxWriter import PdbxWriter
from rcsbapi.model import ModelQuery
from rcsbapi.search import AttributeQuery

class LigandCoordinatesExtract:
    """
    This class provides methods to download and clean mmCIF data for specific ligands in PDB entries.
    It retains only the "chem_comp" and "atom_site" categories.
    """

    def __init__(self):
        self.mQ = ModelQuery()
        self.l_dc = []

    def processOne(self, ccd_id: str, pdb_id: str):
        """
        Retrieve and cleans mmCIF data for a specific ligand in a given PDB entry.
        """
        try:
            # Query the API for atoms related to the ligand in the specified entry
            result = self.mQ.get_atoms(entry_id=pdb_id, label_comp_id=ccd_id)
            if not result:
                print(f"[Warning] No data found for {ccd_id} in {pdb_id}")
                return None

            # Read the result into a file-like object
            file_like = io.StringIO(result)
            pR = PdbxReader(file_like)
            pR.read(self.l_dc, ["chem_comp", "atom_site"])
        except Exception as e:
            print(f"[Error] Failed for {ccd_id} in {pdb_id}: {e}")
            return None
        
    def processAll(self, d_pdb_by_ccd: dict):
        """
        Processes a dictionary mapping ligand IDs to lists of PDB IDs.
        Collects and cleans the relevant mmCIF data for each combination.
        """
        for ccd_id, l_pdb_id in d_pdb_by_ccd.items():
            for pdb_id in l_pdb_id:
                self.processOne(ccd_id, pdb_id)
        
    def writeModel(self, fp_out: str):
        """
        Writes the cleaned mmCIF data to the specified output file.
        """
        try:
            with open(fp_out, "w") as f:
                pW = PdbxWriter(f)
                pW.write(self.l_dc)
            print(f"Output written to {fp_out}")
        except Exception as e:
            print(f"[Error] Failed to write output: {e}")

# def readInputFile(fp_in: str):
#     """
#     Reads the input TSV file and returns a dictionary mapping ligand IDs to lists of PDB IDs.
#     """
#     d_pdb_by_ccd = {}
#     with open(fp_in, "r", encoding="utf-8") as f:
#         for line in f:
#             if line.strip():  # Skip empty lines
#                 ccd_id, pdb_ids = line.strip().split("\t")
#                 d_pdb_by_ccd[ccd_id] = pdb_ids.split()
#     return d_pdb_by_ccd

def searchOneCcd(ccd_id):
    """retrieve all PDB IDs with a specific CCD ID

    Args:
        ccd_id (str): CCD ID, e.g. ATP

    Returns:
        list: PDB IDs list
    """    
    # Build a query for an exact match on the chemical component ID
    query = AttributeQuery(
        attribute="rcsb_chem_comp_container_identifiers.comp_id",  # Path to the comp_id field
        operator="exact_match",                                    # Use exact string matching
        value=ccd_id                                               
    )

    # Execute the query and collect PDB IDs that contain this CCID
    results = list(query(return_type="entry"))  # Get entries (PDB IDs) as the return type
    return [pdb_id.upper() for pdb_id in results]       # Convert results into a set of unique PDB IDs

def searchPdbByCcd(l_ccd, l_pdb_limit=None, n_limit=None):
    """search PDB IDs with the CCD ID

    Args:
        l_ccd (list): CCD IDs to search
        l_pdb_limit (list): limit search to the list of PDB IDs provided
        n_limit (int): limit the number of entries

    Returns:
        d_pdb_by_ccd: dictionary with CCD IDs as keys, and list of PDB IDs as values
    """    
    d_pdb_by_ccd = {}
    for ccd_id in l_ccd:
        d_pdb_by_ccd[ccd_id] = []
        l_pdb = searchOneCcd(ccd_id)
        if l_pdb_limit:
            for pdb_id in l_pdb_limit:
                if pdb_id.upper() in l_pdb:
                    d_pdb_by_ccd[ccd_id].append(pdb_id.upper())
        elif n_limit:
            if len(l_pdb) <= n_limit:
                d_pdb_by_ccd[ccd_id] = l_pdb
            else:
                d_pdb_by_ccd[ccd_id] = l_pdb[:n_limit]
        else:
            d_pdb_by_ccd[ccd_id] = l_pdb

    return d_pdb_by_ccd


def main():
    """
    Retrieve ligand coordinates from PDB entries by author-provided CCD IDs
    """
    parser = argparse.ArgumentParser(description="extract coordinates of specific ligands in PDB entries")

    # take mandatory author-provided CCD IDs, can process multiple IDs
    parser.add_argument("-ccd", dest="ccd", nargs="+", required=True,
                        help="provide ligand CCD ID, e.g. CPT")
    
    # add mutually exclusive options limitating extraction by PDB IDs or number of entries
    limit_extraction = parser.add_mutually_exclusive_group(required=False)
    limit_extraction.add_argument("-pdb", dest="pdb", nargs="+",
                        help="provide PDB IDs to limit the extraction on")
    limit_extraction.add_argument("-n", dest="n",
                        help="run extraction on the limited number of PDB entries")

    args = parser.parse_args()

    l_ccd = args.ccd
    print(f"to extract ligand coordinates for CCD IDs {','.join(l_ccd)}")

    if args.pdb:
        l_pdb = args.pdb
        n = None
        print(f"from PDB IDs {','.join(l_pdb)}")
    elif args.n:
        n = int(args.n)
        l_pdb = None
        print(f"from up to {n} PDB entries")
    else:
        l_pdb = None
        n = None
        print(f"from all PDB entries")

    # fp_in = "cc-to-pdb-picked.tsv"
    # d_pdb_by_ccd = readInputFile(fp_in)

    d_pdb_by_ccd = searchPdbByCcd(l_ccd, l_pdb, n)
  
    if d_pdb_by_ccd:
        for ccd_id, l_pdb in d_pdb_by_ccd.items():
            if l_pdb:
                print(f"to retrieve {ccd_id} coordinates in {len(l_pdb)} PDB entries")
                lCE = LigandCoordinatesExtract()
                for pdb_id in l_pdb:
                    print(f"retrieving {ccd_id} coordinates from {pdb_id}")
                    lCE.processOne(ccd_id, pdb_id)
                fp_out = ccd_id + "-coordinates.cif"
                lCE.writeModel(fp_out)
                print(f"coordinates for {ccd_id} written to {fp_out}")
            else:
                print(f"failed to find PDB entries with CCD ID {ccd_id}")
    else:
        print("Cannot find any PDB entries with provided CCD IDs")

if __name__ == "__main__":
    start = time.time()
    main()
    end = time.time()
    print(f"Processing completed in {end - start:.2f} seconds.")
