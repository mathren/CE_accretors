from MESAreader import getSrcCol
from lib_plot_bin import plot_ratio_BE_r, Rsun_cm
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import glob
import warnings
import paths

rot_root = str(paths.data / "MESA_output/single_stars/Z_0.0019")
grid_folders = sorted(glob.glob(rot_root + "/30_rot0*"))
colors = plt.cm.viridis(np.linspace(0, 1, len(grid_folders)))

fig = plt.figure()
gs = gridspec.GridSpec(100, 100)
ax = fig.add_subplot(gs[:, :])

rot_root = str(paths.data / "MESA_output/single_stars/Z_0.0019")
pfile_nonrot = rot_root + "/30_rot0.0/LOGS/500Rsun.data"
pfiles = sorted(glob.glob(rot_root + "/30_rot0.[123456789]*/LOGS/500Rsun.data"))
colors = plt.cm.viridis(np.linspace(0, 1, len(pfiles)))
# sanity check
ratio = plot_ratio_BE_r(
    pfile_nonrot,
    pfile_nonrot,
    ax,
    alpha_th=1.0,
    alpha_rot=0.0,
    c="orange",
    ls="-",
    lw=2,
    zorder=0,
)
print(min(ratio), max(ratio))
for pfile_rot in pfiles:
    label = pfile_rot.split("/")[-3].split("_rot")[-1]
    label = "$\omega_\mathrm{ZAMS}/\omega_\mathrm{crit}=$" + label
    ratio = plot_ratio_BE_r(
        pfile_rot,
        pfile_nonrot,
        ax,
        alpha_th=1.0,
        alpha_rot=0.0,
        c=colors[pfiles.index(pfile_rot)],
        ls="-",
        lw=2,
        zorder=0,
        label=label
        )
    print(min(ratio), max(ratio))
    # ratio = plot_ratio_BE_r(
    #     pfile_rot,
    #     pfile_nonrot,
    #     ax,
    #     alpha_th=1.0,
    #     alpha_rot=1.0,
    #     c=colors[pfiles.index(pfile_rot)],
    #     ls="--",
    #     lw=5,
    #     zorder=1
    #     )
    # print(min(ratio), max(ratio))
    print("----------------")
ax.axvline(np.log10(500*Rsun_cm), 0,1,ls=":", lw=2, c="#808080", zorder=0)
ax.set_xlim(8.4, 14)
ax.set_ylim(-0.05, 1.75)
ax.set_xlabel(r"$\log_{10}(r/\mathrm{[cm]})$")
ax.set_ylabel(r"BE($\omega$ = 0)/BE($\omega$)")
ax.legend(handlelength=0.5, ncol=2, fontsize=20)
plt.tight_layout()
plt.show()
