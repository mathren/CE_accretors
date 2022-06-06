""" script to test dependence from zenodo data file (intermediate product of processing) """

import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import paths

# try loading
data = np.load(paths.data / "PATH_TO_DATAFILE")

fig = plt.figure()
gs = gridspec.GridSpec(100, 100)
ax = fig.add_subplot(gs[:, :])
ax.plot(np.nan, np.nan, ls="-", c="k", label=r"$^1\mathrm{H}$")
fig.savefig(paths.figures / "try.pdf")
