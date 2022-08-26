from lib_plot_bin import plot_HRD, annotate_radii_hrd
from MESAreader import getSrcCol
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import paths

root = paths.data / "MESA_output"
root_bin = root / "binaries/Z_0.0019"
b1 = (
    root_bin
    / "m1_18.0000_m2_15.0000_initial_z_0.0019_initial_period_in_days_1.0000e+02_grid_index_0_1/"
)


fig = plt.figure(figsize=(20, 20))
gs = gridspec.GridSpec(110, 120)
gs.update(wspace=0, hspace=0)  # top=1.1)
ax = plt.subplot(gs[:, :])

LOGS2 = b1 / "LOGS2/"
hfile2 = str(LOGS2) + "/history.data"



plot_HRD(
    hfile2, ax=ax, annotate_TAMS=False, annotate_RLOF=True, ls="-", c='r', zorder=10
)
# LOGS1 = f / "LOGS1/"
#     hfile1 = str(LOGS1) + "/history.data"
#     plot_HRD(
#         hfile1,
#         ax=ax,
#         annotate_TAMS=False,
#         annotate_RLOF=True,
#         ls="--",
#         lw=2,
#         c=c,
#         zorder=0,
#     )

ax.set_ylim(ymin=4.3, ymax=5.7)
ax.set_xlim(xmin=3.6, xmax=4.75)
annotate_radii_hrd(ax, [10, 50, 100, 200, 300, 500])
ax.invert_xaxis()
# ax.plot(np.nan, np.nan, c="k", ls="--", lw=2, label="donor")
# ax.plot(np.nan, np.nan, c="k", label="accretor")
# ax.legend()
ax.set_xlabel(r"$\log_{10}(T_\mathrm{eff}/[\mathrm{K}])$")
ax.set_ylabel(r"$\log_{10}(L/L_\odot)$")
plt.savefig("/tmp/HRD_single.pdf")
plt.savefig("/tmp/HRD_single.png")
plt.show()
