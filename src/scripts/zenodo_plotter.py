""" script to test dependence from zenodo data file (intermediate product of processing) """

import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import paths
import os

# check if file exists
if os.path.isfile(paths.data / "try.npy"):
    print("try.npy file found")
    print("make dummy fig")

    fig = plt.figure()
    gs = gridspec.GridSpec(100, 100)
    ax = fig.add_subplot(gs[:, :])
    ax.plot(np.linspace(0,1,10), np.linspace(0,1,10), ls="-", c="k")
    fig.savefig(paths.figures / "try.pdf")
