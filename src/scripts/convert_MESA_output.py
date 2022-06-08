"""
Convert the tarball of the entire MESA output to binary.
This allows showyourwork to cache the binary results on zenodo sandbox
for faster reproduction
"""
__author__ = "M. Renzo"
from MESAreader import getSrcCol
import sys
import os
import paths
import glob
import tarfile

if not os.path.isfile(paths.data / "MESA_output.tar.gz"):
    raise FileNotFoundError("MESA output tarball not found! Download from zenodo!")
else:
    tarball = tarfile.open(paths.data / "MESA_output.tar.gz", "r:gz")
    for member in tarball.getmembers():
        # skip folders and only consider *.data files
        if ".data" not in str(member): continue
        # print(member.name)
        bin_fname = paths.data / str(member.name[:-4]+".npy")
        if not os.path.isfile(bin_fname):
            # file does not exist, extract
            tarball.extractall(path = paths.data, members=[member])
            src, col = getSrcCol(paths.data / member.name, True, True, bin_fname=bin_fname)
            # for seeing if all of showyourwork proceeds beyond this script
            os.system("touch "+str(paths.data / "try.npy"))
        else:
            print(str(bin_fname)+" found")
            # for seeing if all of showyourwork proceeds beyond this script
            os.system("touch "+str(paths.data / "try.npy"))
        break
    tarball.close()
