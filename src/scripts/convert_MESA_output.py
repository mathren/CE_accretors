# to be able to split the conversion to *.npy of MESA output for showyourwork
from MESAreader import getSrcCol
import sys
import os
import paths
import glob
import tarfile

if not os.path.isfile(paths.data/"MESA_output.tar.gz"):
    raise FileNotFoundError("MESA output tarball not found! Download from zenodo!")
else:
    tarball = tarfile.open(paths.data/"MESA_output.tar.gz","r:gz")
    for member in tarball.getmembers():
        # skip folders and only consider *.data files
        if ".data" not in str(member):
            pass
        else:
            # print(member)
            f = tarball.extractall(member)
            print(f)
            src, col = getSrcCol(f, True, True)
            break
    tarball.close()



# filenames = glob.glob(paths.data /"MESA_output/**/*.data")
# found = True
# for f in filenames:
#     if not os.path.isfile(f):
#         print("MESA OUTPUT NOT FOUND")
#         found = False
#         break
# if found:
#     os.system("touch src/data/try.npy")
