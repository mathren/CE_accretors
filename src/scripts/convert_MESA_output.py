# to be able to split the conversion to *.npy of MESA output for showyourwork
from MESAreader import getSrcCol
import sys
import os
import paths
import numpy as np
import paths

if __name__ == "__main__":
    input_fname = paths.data / "history_try.data"
    bin_fname = paths.data/"try.npy"
    print("out:", bin_fname)
    if not os.path.isfile(input_fname):
        sys.exit()
    else:
        src, col = getSrcCol(input_fname, True, True, bin_fname)
        np.save(bin_fname, src)
        print("done ", input_fname)
