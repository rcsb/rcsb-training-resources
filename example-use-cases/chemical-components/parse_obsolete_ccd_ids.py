"""
Python script for downloading the current CCD components file and extracting obsolete CCD IDs.
The obsolete IDs are written to a text file, `obsolete_ccd_id.list`.

This script requires the following package, which can be installed with:
pip install mmcif
"""

import os
import gzip
import shutil
import urllib.request
import urllib.error
from mmcif.io.IoAdapterCore import IoAdapterCore


def downloadComponents():
    """Download and unzip the CCD components file from wwPDB.

    Returns:
        str or None: Unzipped local filename on success, otherwise None.
    """
    url = "https://files.wwpdb.org/pub/pdb/data/monomers/components.cif.gz"
    filename_local = url.split('/')[-1]
    print("to download from %s" % url)
    try:
        urllib.request.urlretrieve(url, filename_local)
    except urllib.error.HTTPError as err_http:
        print(err_http)
        print("%s file doesn't exist" % url)
        return None
    except OSError as err_os:
        print(err_os)
        print("cannot write file, check privilege in the current folder")
        return None
    except Exception as err_other:
        print(err_other)
        print("fail to download from %s" % url)
        return None
    else:
        if filename_local.endswith(".gz"):
            print("to unzip %s" % filename_local)
            filename_local_unzipped = filename_local[:-3]  # remove the .gz extension
            try:
                with gzip.open(filename_local, 'rb') as file_in:
                    with open(filename_local_unzipped, 'wb') as file_out:
                        shutil.copyfileobj(file_in, file_out)
                os.remove(filename_local)
                return filename_local_unzipped
            except Exception as err_unzip:
                print(err_unzip)
                print("cannot unzip %s, check file format or privilege in the current folder" % filename_local)
                return None
        else:
            return filename_local

        
def parse(filename_in):
    """Parse a CCD CIF file and collect IDs marked as obsolete.

    Args:
        filename_in (str): Path to the local components CIF file.

    Returns:
        list[str]: List of obsolete CCD IDs.
    """
    l_id = []
    io = IoAdapterCore()
    l_data_container = io.readFile(filename_in)
    for data_container in l_data_container:
        print("check CCD %s" % data_container.getName())
        if data_container.getObj("chem_comp").getValue("pdbx_release_status", 0) == "OBS":
            l_id.append(data_container.getObj("chem_comp").getValue("id", 0))
    return l_id


def write(l_id, filename_out):
    """Write obsolete CCD IDs to a text file, one ID per line.

    Args:
        l_id (list[str]): Obsolete CCD IDs.
        filename_out (str): Output filename.
    """
    with open(filename_out, 'w') as f:
        for ccd_id in l_id:
            f.write(ccd_id)
            f.write("\n")

            
def main():
    """Run the end-to-end workflow for obsolete CCD ID extraction."""
    filename_in = downloadComponents()
    if filename_in is None:
        print("failed to download and/or unzip the components file")
        return
    l_obsolete = parse(filename_in)
    filename_out = "obsolete_ccd_id.list"
    write(l_obsolete, filename_out)
    print("parsed %s obsolete CCD IDs into %s" % (len(l_obsolete), filename_out))


if __name__ == "__main__":
    main()
