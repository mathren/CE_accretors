# to be able to split the conversion to *.npy of MESA output for showyourwork
from MESAreader import getSrcCol
import sys
import os
import paths
import glob

filenames = glob.glob(paths.data /"MESA_output/**/*.data")
for f in filenames:
    print(f)

# input_fname = paths.data / "MESA_output/"




# if not os.path.isfile(input_fname):
#     print("src/data/history.data NOT FOUND")
#     sys.exit()
# else:
#     src, col = getSrcCol(input_fname, True, True, bin_fname)
#     # np.save(bin_fname, src)
#     print("done ", input_fname)
# # if os.path.isfile(paths.data/"MESA_output.tar"):
# #     os.system("touch src/data/try.npy")
# # else:
# #     print("MESA OUTPUT NOT FOUND")
