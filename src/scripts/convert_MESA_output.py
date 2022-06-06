# to be able to split the conversion to *.npy of MESA output for showyourwork
from MESAreader import getSrcCol
import sys
import os
import paths
import numpy as np
import paths

if __name__ == "__main__":
    input_fname = sys.argv[1]
    bin_fname = input_fname.replace('.data', '.npy')
    print("out:", paths.data/bin_fname)

    if not os.path.isfile(input_fname):
        sys.exit()
    else:
        src, col = getSrcCol(input_fname, True, True, paths.data/bin_fname)
        np.save(paths.data / "try.npy", src)
        print("done ", input_fname)
