from MESAreader import getSrcCol
from lib_plot_bin import plot_ratio_BE_r, get_ax_from_pfile, Rsun_cm
from lib_engineered import get_M_boundary
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import glob
import warnings
import paths

def grid_ratios(fig_name=None):
    """
    Parameters:
    ----------
    fig_name: `string` location where to save the figure
    """
    # locations fix with zenodo
    root = paths.data / "MESA_output/"
    root_eng = root / "engineered_stars/same_core/"
    grid_folders18 = sorted(glob.glob(str(root_eng) + "/grid18/*.*/LOGS/"))
    grid_folders20 = sorted(glob.glob(str(root_eng) + "/grid20/*.*/LOGS/"))
    grid_folders36 = sorted(glob.glob(str(root_eng) + "/grid36/*.*/LOGS/"))
    ## accretor models
    root_accretors = root / "binaries/Z_0.0019/"
    #
    b1 = (
        str(root_accretors)
        + "/m1_18.0000_m2_15.0000_initial_z_0.0019_initial_period_in_days_1.0000e+02_grid_index_0_1/LOGS2/"
    )
    b2 = (
        str(root_accretors)
        + "/m1_20.0000_m2_17.0000_initial_z_0.0019_initial_period_in_days_1.0000e+02_grid_index_0_1/LOGS2/"
    )
    b3 = (
        str(root_accretors)
        + "/m1_38.0000_m2_30.0000_initial_z_0.0019_initial_period_in_days_1.0000e+02_grid_index_0_1/LOGS2/"
    )

    single_grid1 = grid_folders18
    single_grid2 = grid_folders20
    single_grid3 = grid_folders36

    s1 = str(root) + "/single_stars/Z_0.0019/18_rot0.0/LOGS/"
    s2 = str(root) + "/single_stars/Z_0.0019/20_rot0.0/LOGS/"
    s3 = str(root) + "/single_stars/Z_0.0019/36_rot0.0/LOGS/"

    fig = plt.figure(figsize=(16, 24))
    gs = gridspec.GridSpec(120, 150)
    # for sanity checks
    fig2, dummy_ax = plt.subplots()
    dummy_ax.set_ylim(-0.1, 0.1)
    # accretor15, compare to 18Msun
    ax1 = fig.add_subplot(gs[:20, :50])
    ax2 = fig.add_subplot(gs[20:40, :50])
    ax3 = fig.add_subplot(gs[40:60, :50])
    ax4 = fig.add_subplot(gs[60:80, :50])
    ax5 = fig.add_subplot(gs[80:100, :50])

    axes = [ax1, ax2, ax3, ax4, ax5]
    for ax in axes:
        ax.set_xlim(8.2, 14)
        ax.set_ylim(-0.05, 2.5)
        if ax != axes[-1]:
            ax.set_xticklabels([])
            ax.axhline(1, 0, 1, lw=2, ls="--", c="#808080")
            ax.text(
                0.05,
                0.85,
                "accretor more bound",
                fontsize=20,
                transform=ax.transAxes,
                va="bottom",
                ha="left",
                zorder=0,
            )
            ax.text(
                0.05,
                0.05,
                "accretor less bound",
                fontsize=20,
                transform=ax.transAxes,
                va="bottom",
                ha="left",
                zorder=0,
            )
        else:
            ax.set_xlabel(r"$\log_{10}(r/\mathrm{[cm]})$")
    ax3.set_ylabel(r"$BE(\mathrm{accretor})/BE(\mathrm{single})$")

    profiles = sorted(glob.glob(b1 + "/*Rsun.data"))
    for pfile_accretor in profiles:
        # # sanity checks
        # ratio = plot_ratio_BE_r(pfile_accretor, pfile_accretor, dummy_ax, alpha_th=0)
        # if (round(max(ratio), 5) != 1 or round(min(ratio), 5) != 1):
        #     print("0", pfile_accretor, np.nanmax(ratio), np.nanmin(ratio))
        # ratio = plot_ratio_BE_r(pfile_accretor, pfile_accretor, dummy_ax, alpha_th=1)
        # if (round(max(ratio), 5) != 1 or round(min(ratio), 5) != 1):
        #     print("1", pfile_accretor, np.nanmax(ratio), np.nanmin(ratio))
        # # -----------------------
        string = pfile_accretor.split("/")[-1]
        ax = get_ax_from_pfile(pfile_accretor, axes)
        colors = plt.cm.viridis(np.linspace(0, 1, len(single_grid1)))
        for f in single_grid1:
            pfile_single = glob.glob(f + string)[0]  # +string)
            # print(colored(pfile_single, "green"))
            # gravitational only
            # alpha_th = 0.0
            # # sanity check
            # # ratio = plot_ratio_BE_r(pfile_single, pfile_single, dummy_ax, alpha_th=alpha_th)
            # # if (round(max(ratio), 5) != 1 or round(min(ratio), 5) != 1):
            # #     print("0", pfile_single, np.nanmax(ratio), np.nanmin(ratio))
            # # accretor/single
            # #################
            # ratio = plot_ratio_BE_r(
            #     pfile_single,
            #     pfile_accretor,
            #     ax,
            #     alpha_th=alpha_th,
            #     color=colors[single_grid1.index(f)],
            # )
            #################
            # internal energy included
            alpha_th = 1.0
            # ratio = plot_ratio_BE_r(pfile_single, pfile_single, dummy_ax, alpha_th=alpha_th)
            # if (round(max(ratio), 5) != 1 or round(min(ratio), 5) != 1):
            #     print("0", pfile_single, np.nanmax(ratio), np.nanmin(ratio))
            # # accretor/single
            ratio = plot_ratio_BE_r(
                pfile_single,
                pfile_accretor,
                ax,
                alpha_th=alpha_th,
                color=colors[single_grid1.index(f)],
                ls="-", lw=2
            )
        # now plot normal single star
        pfile_single = glob.glob(s1 + string)[0]
        # alpha_th = 0.0
        # # sanity check
        # # ratio = plot_ratio_BE_r(pfile_single, pfile_single, dummy_ax, alpha_th=alpha_th)
        # # if (round(max(ratio), 5) != 1 or round(min(ratio), 5) != 1):
        # #     print("0", pfile_single, np.nanmax(ratio), np.nanmin(ratio))
        # # accretor/single
        # #################
        # ratio = plot_ratio_BE_r(
        #     pfile_single, pfile_accretor, ax, alpha_th=alpha_th, color="r"
        # )
        # #################
        # internal energy included
        alpha_th = 1.0
        # ratio = plot_ratio_BE_r(pfile_single, pfile_single, dummy_ax, alpha_th=alpha_th)
        # if (round(max(ratio), 5) != 1 or round(min(ratio), 5) != 1):
        #     print("0", pfile_single, np.nanmax(ratio), np.nanmin(ratio))
        # # accretor/single
        ratio = plot_ratio_BE_r(
            pfile_single, pfile_accretor,  ax, alpha_th=alpha_th, color="r", ls="-"
        )
        # print(pfile_single, min(ratio), max(ratio))
    # # 20 Msun
    bx1 = fig.add_subplot(gs[:20, 50:100])
    bx2 = fig.add_subplot(gs[20:40, 50:100])
    bx3 = fig.add_subplot(gs[40:60, 50:100])
    bx4 = fig.add_subplot(gs[60:80, 50:100])
    bx5 = fig.add_subplot(gs[80:100, 50:100])

    bxes = [bx1, bx2, bx3, bx4, bx5]
    for bx in bxes:
        bx.set_xlim(8.2, 14)
        bx.set_ylim(-0.05, 2.5)
        bx.set_yticklabels([])
        if bx != bxes[-1]:
            bx.set_xticklabels([])
            bx.axhline(1, 0, 1, lw=2, ls="--", c="#808080")
            bx.text(
                0.05,
                0.85,
                "accretor more bound",
                fontsize=20,
                transform=bx.transAxes,
                va="bottom",
                ha="left",
                zorder=0,
            )
            bx.text(
                0.05,
                0.05,
                "accretor less bound",
                fontsize=20,
                transform=bx.transAxes,
                va="bottom",
                ha="left",
                zorder=0,
            )
        else:
            bx.set_yticklabels([])
            bx.set_xlabel(r"$\log_{10}(r/\mathrm{[cm]})$")

    profiles = sorted(glob.glob(b2 + "/*Rsun.data"))  # , key=get_age_from_profile)
    # for pfile_accretor in tqdm(profiles, desc="20Msun"):
    for pfile_accretor in profiles:
        # print()
        # print(colored(pfile_accretor, "blue"))
        # # sanity checks
        # ratio = plot_ratio_BE_r(pfile_accretor, pfile_accretor, dummy_ax, alpha_th=0)
        # if (round(max(ratio), 5) != 1 or round(min(ratio), 5) != 1):
        #     print("0", pfile_accretor, np.nanmax(ratio), np.nanmin(ratio))
        # ratio = plot_ratio_BE_r(pfile_accretor, pfile_accretor, dummy_ax, alpha_th=1)
        # if (round(max(ratio), 5) != 1 or round(min(ratio), 5) != 1):
        #     print("1", pfile_accretor, np.nanmax(ratio), np.nanmin(ratio))
        # -----------------------
        string = pfile_accretor.split("/")[-1]
        bx = get_ax_from_pfile(pfile_accretor, bxes)
        colors = plt.cm.viridis(np.linspace(0, 1, len(single_grid2)))
        for f in single_grid2:
            pfile_single = glob.glob(f + string)[0]  # +string)
            # print(colored(pfile_single, "green"))
            # gravitational only
            # alpha_th = 0.0
            # # sanity check
            # # ratio = plot_ratio_BE_r(pfile_single, pfile_single, dummy_ax, alpha_th=alpha_th)
            # # if (round(max(ratio), 5) != 1 or round(min(ratio), 5) != 1):
            # #     print("0", pfile_single, np.nanmax(ratio), np.nanmin(ratio))
            # # accretor/single
            # #################
            # ratio = plot_ratio_BE_r(
            #     pfile_single,
            #     pfile_accretor,
            #     bx,
            #     alpha_th=alpha_th,
            #     color=colors[single_grid2.index(f)],
            # )
            #################
            # internal energy included
            alpha_th = 1.0
            # ratio = plot_ratio_BE_r(pfile_single, pfile_single, dummy_ax, alpha_th=alpha_th)
            # if (round(max(ratio), 5) != 1 or round(min(ratio), 5) != 1):
            #     print("0", pfile_single, np.nanmax(ratio), np.nanmin(ratio))
            # # accretor/single
            ratio = plot_ratio_BE_r(
                pfile_single,
                pfile_accretor,
                bx,
                alpha_th=alpha_th,
                color=colors[single_grid2.index(f)],
                ls="-", lw=2
            )
        # now plot normal single star
        pfile_single = glob.glob(s2 + string)[0]
        # alpha_th = 0.0
        # # sanity check
        # # ratio = plot_ratio_BE_r(pfile_single, pfile_single, dummy_ax, alpha_th=alpha_th)
        # # if (round(max(ratio), 5) != 1 or round(min(ratio), 5) != 1):
        # #     print("0", pfile_single, np.nanmax(ratio), np.nanmin(ratio))
        # # accretor/single
        # #################
        # ratio = plot_ratio_BE_r(
        #     pfile_single, pfile_accretor, bx, alpha_th=alpha_th, color="r"
        # )
        #################
        # internal energy included
        alpha_th = 1.0
        # ratio = plot_ratio_BE_r(pfile_single, pfile_single, dummy_ax, alpha_th=alpha_th)
        # if (round(max(ratio), 5) != 1 or round(min(ratio), 5) != 1):
        #     print("0", pfile_single, np.nanmax(ratio), np.nanmin(ratio))
        # # accretor/single
        ratio = plot_ratio_BE_r(
            pfile_single,
            pfile_accretor,
            bx,
            alpha_th=alpha_th,
            color="r",
            ls="-",
        )
        print(pfile_single, min(ratio), max(ratio))
    # 30
    cx1 = fig.add_subplot(gs[:20, 100:])
    cx2 = fig.add_subplot(gs[20:40, 100:])
    cx3 = fig.add_subplot(gs[40:60, 100:])
    cx4 = fig.add_subplot(gs[60:80, 100:])
    cx5 = fig.add_subplot(gs[80:100, 100:])

    cxes = [cx1, cx2, cx3, cx4, cx5]
    for cx in cxes:
        cx.set_xlim(8.2, 14)
        cx.set_ylim(-0.05, 2.5)
        cx.text(
            0.05,
            0.85,
            "accretor more bound",
            fontsize=20,
            transform=cx.transAxes,
            va="bottom",
            ha="left",
            zorder=0,
        )
        cx.text(
            0.05,
            0.05,
            "accretor less bound",
            fontsize=20,
            transform=cx.transAxes,
            va="bottom",
            ha="left",
            zorder=0,
        )
        cx.axhline(1, 0, 1, lw=2, ls="--", c="#808080")
        cx.set_yticklabels([])
        if cx != cxes[-1]:
            cx.set_xticklabels([])
        else:
            cx.set_xlabel(r"$\log_{10}(r/\mathrm{[cm]})$")

    profiles = sorted(glob.glob(b3 + "/*Rsun.data"))  # , key=get_age_from_profile)
    # for pfile_accretor in tqdm(profiles, desc="30Msun"):
    for pfile_accretor in profiles:
        # print()
        # print(colored(pfile_accretor, "blue"))
        # # sanity checks
        # ratio = plot_ratio_BE_r(pfile_accretor, pfile_accretor, dummy_ax, alpha_th=0)
        # if (round(max(ratio), 5) != 1 or round(min(ratio), 5) != 1):
        #     print("0", pfile_accretor, np.nanmax(ratio), np.nanmin(ratio))
        # ratio = plot_ratio_BE_r(pfile_accretor, pfile_accretor, dummy_ax, alpha_th=1)
        # if (round(max(ratio), 5) != 1 or round(min(ratio), 5) != 1):
        #     print("1", pfile_accretor, np.nanmax(ratio), np.nanmin(ratio))
        # -----------------------
        string = pfile_accretor.split("/")[-1]
        cx = get_ax_from_pfile(pfile_accretor, cxes)
        colors = plt.cm.viridis(np.linspace(0, 1, len(single_grid3)))
        for f in single_grid3:
            # print(colored(pfile_single, "green"))
            pfile_single = glob.glob(f + string)[0]  # +string)
            # gravitational only
            # alpha_th = 0.0
            # # sanity check
            # # ratio = plot_ratio_BE_r(pfile_single, pfile_single, dummy_ax, alpha_th=alpha_th)
            # # if (round(max(ratio), 5) != 1 or round(min(ratio), 5) != 1):
            # #     print("0", pfile_single, np.nanmax(ratio), np.nanmin(ratio))
            # # accretor/single
            # #################
            # ratio = plot_ratio_BE_r(
            #     pfile_single,
            #     pfile_accretor,
            #     cx,
            #     alpha_th=alpha_th,
            #     color=colors[single_grid3.index(f)],
            # )
            #################
            # internal energy included
            alpha_th = 1.0
            # ratio = plot_ratio_BE_r(pfile_single, pfile_single, dummy_ax, alpha_th=alpha_th)
            # if (round(max(ratio), 5) != 1 or round(min(ratio), 5) != 1):
            #     print("0", pfile_single, np.nanmax(ratio), np.nanmin(ratio))
            # accretor/single
            ratio = plot_ratio_BE_r(
                pfile_single,
                pfile_accretor,
                cx,
                alpha_th=alpha_th,
                color=colors[single_grid3.index(f)],
                ls="-", lw=2
            )
        # now plot normal single star
        pfile_single = glob.glob(s3 + string)[0]
        # alpha_th = 0.0
        # # sanity check
        # # ratio = plot_ratio_BE_r(pfile_single, pfile_single, dummy_ax, alpha_th=alpha_th)
        # # if (round(max(ratio), 5) != 1 or round(min(ratio), 5) != 1):
        # #     print("0", pfile_single, np.nanmax(ratio), np.nanmin(ratio))
        # # accretor/single
        # #################
        # ratio = plot_ratio_BE_r(
        #     pfile_single, pfile_accretor, cx, alpha_th=alpha_th, color="r"
        # )
        # #################
        # internal energy included
        alpha_th = 1.0
        # ratio = plot_ratio_BE_r(pfile_single, pfile_single, dummy_ax, alpha_th=alpha_th)
        # if (round(max(ratio), 5) != 1 or round(min(ratio), 5) != 1):
        #     print("0", pfile_single, np.nanmax(ratio), np.nanmin(ratio))
        # # accretor/single
        ratio = plot_ratio_BE_r(
            pfile_single, pfile_accretor, cx, alpha_th=alpha_th, color="r", ls="-"
        )
        print(pfile_single, min(ratio), max(ratio))
    for cx in cxes:
        dx = cx.twinx()
        dx.set_yticks(cx.get_yticks())
        dx.set_yticklabels([])
        if cx == cxes[0]:
            # bx.text(0.2, 0.1, , fontsize=30, transform=bx.transAxes, va="bottom", ha="right")
            label = "$100\,R_\odot$"
            X = 100 * Rsun_cm
        if cx == cxes[1]:
            # bx.text(0.2, 0.1, , fontsize=30, transform=bx.transAxes, va="bottom", ha="right")
            label = "$200\,R_\odot$"
            X = 200 * Rsun_cm
        if cx == cxes[2]:
            # bx.text(0.2, 0.1, , fontsize=30, transform=bx.transAxes, va="bottom", ha="right")
            label = "$300\,R_\odot$"
            X = 300 * Rsun_cm
        if cx == cxes[3]:
            # bx.text(0.2, 0.1, , fontsize=30, transform=bx.transAxes, va="bottom", ha="right")
            label = "$500\,R_\odot$"
            X = 500 * Rsun_cm
        if cx == cxes[4]:
            # bx.text(0.2, 0.1, , fontsize=30, transform=bx.transAxes, va="bottom", ha="right")
            label = "$1000\,R_\odot$"
            X = 1000 * Rsun_cm
        dx.set_ylabel(label)
        cx.axvline(np.log10(X), 0, 1, ls="--", lw=2, c="#808080", zorder=0)
        # bx.axhline(0, 0,1, ls='--', lw=2, c="#808080", zorder=0)
        ax = axes[cxes.index(cx)]
        if ax != axes[-1]:  # skip last panel
            ax.axvline(np.log10(X), 0, 1, ls="--", lw=2, c="#808080", zorder=0)
            # ax.axhline(0, 0,1, ls='--', lw=2, c="#808080", zorder=0)
        bx = bxes[cxes.index(cx)]
        if bx != bxes[-1]:  # skip last panel
            bx.axvline(np.log10(X), 0, 1, ls="--", lw=2, c="#808080", zorder=0)

    ax1.set_title(r"$M_2=15\rightarrow 18\,M_\odot$", size=30)
    bx1.set_title(r"$M_2=17\rightarrow 20\,M_\odot$", size=30)
    cx1.set_title(r"$M_2=30\rightarrow 36\,M_\odot$", size=30)

    # legends
    # bxes[-1].plot(np.nan, np.nan, c="k", ls="--", label=r"$\alpha_\mathrm{th}=0$")
    # bxes[-1].plot(np.nan, np.nan, c="k", ls="-", label=r"$\alpha_\mathrm{th}=1$")
    # bxes[-1].legend(loc="center")

    axes[-1].plot(np.nan, np.nan, c="r", ls="-", label="ratio to single")
    axes[-1].legend(handletextpad=0.1, handlelength=0.75, loc="center")

    plt.close(fig2)
    if fig_name:
        plt.savefig(fig_name)


if __name__ == "__main__":
    # three_masses_grid(fig_name='grid_ratios.pdf')
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        grid_ratios(fig_name=paths.figures / "grid_ratios.pdf")
