# Author:
#          Mathieu Renzo <mrenzo@flatironinstitute.org>
#
# Keywords: files

# Copyright (C) 2020-2021 Mathieu Renzo

# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or (at
# your option) any later version.
#
# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see http://www.gnu.org/licenses/.

import sys
import socket
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import numpy as np
import re
import os
import warnings
import glob

try:
    from termcolor import colored
except ImportError:

    def colored(a, color):
        return a


from lib_engineered import get_M_boundary

try:
    from MESAreader import getSrcCol, secyer, G_cgs, Lsun, Msun, Rsun_cm, clight
except:
    raise ModuleNotFoundError("MESAreader")

from scipy.interpolate import interp1d


# -------------------------------------------------------------------
# HRD
def plot_HRD(input_file, ax, convert=True, annotate_radii=None, **plot_kwargs):
    """Make HRD

    Parameters:
    ----------
    input_file: `string` or `os.path` of the history.data file
    convert: `bool` optional, create *.npy binary of the data file
    ax: `mpl.axes`
    annotate_radii: `None` or `np.array(dtype=float)`  optional, annotate loci of constant R provided
    """
    # print(colored(input_file, "blue"))
    if os.path.isdir(input_file):
        src, col = getSrcCol(input_file + "/history.data", convert, convert)
    elif os.path.isfile(input_file):
        src, col = getSrcCol(input_file, convert, convert)
    else:
        raise FileNotFoundError(input_file)
    logT = src[:, col.index("log_Teff")]
    logL = src[:, col.index("log_L")]
    try:
        # annotate masses
        fff = input_file[:-7]
        m1, m2 = get_masses(fff)
        if "LOGS1" in input_file:
            ax.text(
                logT[0] - 0.01,
                logL[0] - 0.01,
                f"{m1:.1f}",
                fontsize=20,
                va="center",
                ha="center",
                transform=ax.transData,
            )
        elif "LOGS2" in folder:
            ax.text(
                logT[0] - 0.01,
                logL[0] - 0.01,
                f"{m2:.1f}",
                fontsize=20,
                va="center",
                ha="center",
                transform=ax.transData,
            )
    except:
        pass
    ax.plot(logT, logL, **plot_kwargs)
    if annotate_radii:
        annotate_radii_hrd(ax, radii=annotate_radii)


def get_L_from_r_teff(radius, teff):
    # to annotate radii on HRD
    from math import pi

    # Stephan Boltzman constant
    boltzm = 1.380649e-16  # cgs
    hbar = 6.62607015e-27 / (2 * pi)
    clight = 2.99792458e10
    sigma = (pi * pi * boltzm * boltzm * boltzm * boltzm) / (
        60 * hbar * hbar * hbar * clight * clight
    )
    # convert r to cm
    radius *= Rsun_cm
    # assume teff is in K
    l = 4 * pi * radius * radius * sigma * teff ** 4.0
    # convert to Lsun
    l = l / Lsun
    return l


def annotate_radii_hrd(ax, radii=np.logspace(0, 3, base=10)):
    """
    give the axis object for an HRD plot (assumed to be in log10 Lsun and log10 Teff),
    and a list of radii in Rsun units, plots radii.
    Parameters:
    ----------
    ax: `mpl.ax` matplotlib axis object
    radii: `np.array`, optional, radii to mark
    """
    xmin, xmax = ax.get_xlim()
    ymin, ymax = ax.get_ylim()
    Tmax = 10.0 ** xmax
    Tmin = 10.0 ** xmin
    teff = np.linspace(Tmin, Tmax, 10)
    x = np.log10(teff)
    for r in radii:
        l = get_L_from_r_teff(r, teff)
        y = np.log10(l)
        ax.plot(x, y, c="#808080", ls="-.", lw=1)
        # ax.text(x[5], y[5], f"{r:.0f}"+r"$\,R_\odot$", fontsize=20, transform=ax.transData, zorder=0, c="#808080")# , rotation=np.((max(y)-min(y))/(max(x)-min(x))))
    # reset ylim
    ax.set_ylim(ymin, ymax)


# -------------------------------------------------------------------
# orbital velocity evolution
# def plot_orbital_v2(ax, bfile, c="#77CCCC", label=""):
#     srcb, colb = getSrcCol(bfile)
#     v2 = srcb[:, colb.index("v_orb_2")]
#     t = srcb[:, colb.index("age")] * 1e-6
#     ax.plot(t, v2, lw=3, c=c, label=label)
#     ax.text(t[-1] + 0.1, v2[-1], f"{v2[-1]:.1f}", fontsize=30)


# -------------------------------------------------------------------
# chemical composition
# def get_epsilon_he(X, Y, Mtot=20, dq=1e-8):
#     """
#     convert mass fractions from MESA to the epslion spectroscopists like

#     inputs:

#     X, Y mass fractions of hydrogen and helium respectively,
#     Mtot total mass
#     dq fraction of the total mass that is the surface, 1e-8 is the MESA default

#     output: epsilon_he = N(He)/(N(H)+N(He)) with N number abundance
#     """
#     # From $MESA_DIR/const/public/const_def.f90
#     mp = 1.67262192369e-24  # grams
#     mhe = 4 * mp  # neglects the binding energy of the he4 atom
#     # get mass of the photosphere in grams
#     mass = Mtot * Msun * dq
#     # get number abundances from mass and mass fractions
#     N_h = mass * X / mp
#     N_he = mass * Y / mhe
#     # calculate epsilon
#     epsilon = N_he / (N_h + N_he)
#     return epsilon


# def get_epsilon(mass_frac, X):
#     """
#     convert mass fractions from MESA to the epslion spectroscopists like

#     inputs:

#     X mass fraction of hydrogen
#     mass_frac mass fractions of the wanted element

#     output: 12+log(mass_frac/X)

#     """
#     return 12 + np.log10(mass_frac / X)


def plot_surface_abundances(hfile1, ax=None, label=None, legend=False, fig_name=None):
    """
    plot the surface abundances of a few isotopes
    the post binary evolution is optional

    dash-dotted lines are the initial abundances.
    solid lines are the current abundances.
    dashed line are the conversion of the epsilon
    values reported in the literature for Zeta Ophiuchi. These depend
    on the amount of Y, so they have the same time dependence of the
    abundance of helium. Data from Villamariz & Herrero 2005.
    """
    src, col = getSrcCol(hfile1)

    he4_1 = src[:, col.index("surface_he4")]
    h1_1 = src[:, col.index("surface_h1")]
    n14_1 = src[:, col.index("surface_n14")]
    c12_1 = src[:, col.index("surface_c12")]
    o16_1 = src[:, col.index("surface_o16")]
    t_1 = src[:, col.index("star_age")] * 1e-6

    he4 = he4_1
    h1 = h1_1
    n14 = n14_1
    c12 = c12_1
    o16 = o16_1
    t = t_1

    if ax == None:
        # create plot figure if ax is not passed
        fig = plt.figure()
        gs = gridspec.GridSpec(100, 100)
        ax = fig.add_subplot(gs[:, :])

    if label != None:
        ax.set_title(label, fontsize=30)

    ax.plot(t, h1, c="b", label=r"$^1\mathrm{H}$")
    ax.plot(t, he4, c="r", label=r"$^4\mathrm{He}$")
    ax.plot(t, c12, c="g", label=r"$^{12}\mathrm{C}$")
    ax.plot(t, n14, c="m", label=r"$^{14}\mathrm{N}$")
    ax.plot(t, o16, c="y", label=r"$^{16}\mathrm{O}$")

    ax.axhline(h1[0], 0, 1, c="b", ls="-.", lw=1)
    ax.axhline(he4[0], 0, 1, c="r", ls="-.", lw=1)
    ax.axhline(n14[0], 0, 1, c="m", ls="-.", lw=1)
    ax.axhline(c12[0], 0, 1, c="g", ls="-.", lw=1)
    ax.axhline(o16[0], 0, 1, c="y", ls="-.", lw=1)

    ax.set_xlabel(r"$\mathrm{time \ [Myr]}$")
    ax.set_ylabel(r"$\mathrm{Surface\ mass\ fraction}\ X_i$")
    ax.set_yscale("log")
    if legend:
        ax.legend(ncol=2, handletextpad=0.5, handlelength=0.5, columnspacing=0.75)
    if fig_name:
        plt.savefig(fig_name)
    # return the mass fractions at the end of the run
    return (h1[-1], he4[-1], n14[-1], c12[-1], o16[-1])


# def MassVelocityEvolution(folder, convert=False, figName=""):
#     """folder is the binary evolution MESA folder"""
#     import matplotlib.pyplot as plt
#     import matplotlib.gridspec as gridspec

#     print(folder)
#     fig = plt.figure(figsize=(15, 9))
#     gs = gridspec.GridSpec(200, 100)
#     ax = fig.add_subplot(gs[:100, :])
#     bx = fig.add_subplot(gs[100:200, :])
#     bbx = bx.twinx()

#     srcb, colb = getSrcCol(folder + "/binary_history.data", convert, convert)
#     t = srcb[:, colb.index("age")] * 1e-6
#     M1 = srcb[:, colb.index("star_1_mass")]
#     M2 = srcb[:, colb.index("star_2_mass")]
#     v2 = srcb[:, colb.index("v_orb_2")]
#     P = srcb[:, colb.index("period_days")]

#     ax.plot(t, M2, c="r", label=r"$M_2$")
#     ax.plot(t, M1, c="b", label=r"$M_1$")
#     ax.plot(t, M1 + M2, c="k", label=r"$M_1+M_2$")

#     ax.text(t[-1] + 0.1, M2[-1], r"$M_2=" + f"{M2[-1]:.1f}" + r"$", fontsize=30)
#     ax.text(t[-1] + 0.1, M1[-1], r"$M_1=" + f"{M1[-1]:.1f}" + r"$", fontsize=30)
#     ax.text(t[-1] + 0.1, M1[-1] + M2[-1], r"$M_1+M_2=" + f"{M1[-1]+M2[-1]:.1f}" + r"$", fontsize=30)

#     bx.plot(t, v2, ls="-", lw=3, c="r")
#     bx.text(t[-1] + 0.1, v2[-1], f"{v2[-1]:.1f}", fontsize=30)
#     bbx.plot(t, P, ls="--", lw=3, c="b")

#     # ax.set_ylim(0,37)
#     ax.set_xlim(0, 12)
#     ax.set_xticklabels([])
#     bx.set_xlim(ax.get_xlim())
#     # bbx.set_ylim(90,700)
#     bx.set_xlabel(r"$\mathrm{t \ [Myr]}$")
#     ax.set_ylabel(r"$M \ [M_\odot]$")
#     bx.set_ylabel(r"$v_2 \ [\mathrm{km\ s^{-1}}]$")
#     bbx.set_ylabel(r"$P \ \mathrm{[days]}$", color="b")

#     if figName != "":
#         plt.savefig(figName)


# def get_age_from_profile(pfile):
#     # print(pfile)
#     with open(pfile, 'r') as f:
#         for i, line in enumerate(f):
#             if i == 1:
#                 header_cols = line.split()
#                 # print(header_cols)
#             if i == 2:
#                 header_data = line.split()
#                 break
#     age = float(header_data[header_cols.index('star_age')]) * 1e-6
#     return age  # in Myr


# def get_modnum_from_profile(pfile):
#     # print(pfile)
#     with open(pfile, 'r') as f:
#         for i, line in enumerate(f):
#             if i == 1:
#                 header_cols = line.split()
#                 # print(header_cols)
#             if i == 2:
#                 header_data = line.split()
#                 break
#     mn = int(header_data[header_cols.index('model_number')])
#     return mn  # in yr


# def get_ZAMS_abundances(hfile):
#     src, col = getSrcCol(hfile)
#     surface_c12 = src[0, col.index("surface_c12")]
#     surface_n14 = src[0, col.index("surface_n14")]
#     surface_o16 = src[0, col.index("surface_o16")]
#     return surface_c12, surface_n14, surface_o16


def get_masses(f):
    """returns the initial masses of a MESA run reading them from the folder name
    It assumes the folder names follow the POSYDON convention

    f: string with the absolute path of the working directory
    output: m1, m2 (or m1, nan for a single star) as floats
    """
    if f[-1] == "/":
        folder = f.split("/")[-2]
    elif f[-1] != "/":
        folder = f.split("/")[-1]
    # print(folder)
    try:
        m = re.findall("[+-]?\d+\.\d+", folder)
        # print(m)
        if len(m) == 1:
            # print("single star!")
            m1 = float(m[0])
            m2 = np.nan
        else:
            # print("binary!")
            m1 = float(m[0])
            m2 = float(m[1])
    except:
        # read from src assuming single star
        src, col = getSrcCol(f + "LOGS/history.data", False, False)
        m1 = src[0, col.index("star_mass")]
        ms = np.nan
    return m1, m2


# def plot_rho_mass(pfile, ax="", **plot_kwargs):
#     print(colored(pfile, "blue"))
#     if ax == "":
#         fig = plt.figure()
#         gs = gridspec.GridSpec(100, 110)
#         ax = fig.add_subplot(gs[:, :])
#     src, col = getSrcCol(pfile)
#     logRho = src[:, col.index("logRho")]
#     m = src[:, col.index("mass")]
#     ax.plot(m, logRho, **plot_kwargs)


# def plot_mdot_t(hfile, ax="", **plot_kwargs):
#     print(colored(hfile, "blue"))
#     src, col = getSrcCol(hfile)
#     log_abs_mdot = src[:, col.index("log_abs_mdot")]
#     t = src[:, col.index("star_age")] * 1e-6
#     if ax == "":
#         fig = plt.figure()
#         gs = gridspec.GridSpec(100, 110)
#         ax = fig.add_subplot(gs[:, :])
#     ax.plot(t, log_abs_mdot, **plot_kwargs)


# def plot_omega_mass(pfile, ax="", bx="", **plot_kwargs):
#     print(colored(pfile, "blue"))
#     if ax == "":
#         fig = plt.figure()
#         gs = gridspec.GridSpec(100, 110)
#         ax = fig.add_subplot(gs[:, :])
#     src, col = getSrcCol(pfile)
#     omega = src[:, col.index("omega")] * 60 * 60 * 24
#     m = src[:, col.index("mass")]
#     ax.plot(m, omega, c=c, **plot_kwargs)
#     if bx != "":
#         c12 = src[:, col.index("c12")]
#         bx.plot(m, c12, **plot_kwargs)


# def get_t_free_fall(pfile):
#     """dynamic timescale (seconds) -- estimated by 2*pi*sqrt(r^3/(G*m))"""
#     src, col = getSrcCol(pfile)
#     radius = src[:, col.index("radius")] * Rsun_cm  # cm
#     mass = src[:, col.index("mass")] * Msun  # in g
#     t_ff = 2 * math.pi * np.sqrt(radius * radius * radius / (G_cgs * mass))
#     return t_ff  # in sec


# def get_omega_div_omega_crit(folder):
#     """
#     get the initial rotation rate from the history file
#     The history file is assumed to be in folder/LOGS/history.data
#     """
#     hfile = folder+'LOGS/history.data'
#     src, col = getSrcCol(hfile)
#     omega_div_omega_crit = src[0, col.index("surf_avg_omega_div_omega_crit")]
#     return round(omega_div_omega_crit,1)

# def plot_moment_of_inertia(pfile, **plot_kwargs):
#     print(colored(pfile, "blue"))
#     src, col = getSrcCol(pfile)
#     m = src[:, col.index("mass")]
#     try:
#         i = src[:, col.index("i_rot")]
#     except:
#         j = src[:, col.index("j_rot")]
#         omega = src[:, col.index("omega")]
#         i = j / omega
#     ax.plot(m, i, **plot_kwargs)


# def get_Z_from_hfile(hfile):
#     """
#     checks the initial metallicity from surface at beginning of
#     history.data
#     """
#     src, col = getSrcCol(hfile)
#     # get surface hydrogen
#     X = src[0, col.index('surface_h1')]
#     try:
#         X += src[0, col.index('surface_prot')]
#     except:
#         pass
#         # print(colored("no prot","yellow"))
#     try:
#         X += src[0, col.index('surface_h2')]
#     except:
#         pass
#         # print(colored("no h2","yellow"))
#     # get surface helium
#     Y = src[0, col.index('surface_he4')]
#     try:
#         Y += src[0, col.index('surface_he3')]
#     except:
#         pass
#         # print(colored("no he3","yellow"))
#     Z = 1.0 - X - Y
#     # print(f"Initial Z:{Z:.4f}")
#     return Z


# def get_folder_from_file(fname):
#     """given the filename get profile or history of parent folder"""
#     folder_name = fname.split('/')[-3]
#     return folder_name


# def get_resolution_from_file(fname):
#     """given path of pfile or hfile get resolution"""
#     folder_name = get_folder_from_file(fname)
#     resolution = "$" + (folder_name.split('e_')[-1].split('_J')[0]).replace('_', '') + "$"
#     return resolution


# def plot_n_modnum(hfile, ax="", bx="", **plot_kwargs):
#     src, col = getSrcCol(hfile)
#     nz = src[:, col.index("num_zones")]
#     modnum = src[:, col.index("model_number")]
#     log_dt = src[:, col.index("log_dt")]
#     tc = src[:, col.index("log_center_T")]
#     i = np.argmin(np.absolute(tc-8.2))
#     ax.axvline(modnum[i], 0,1,lw=1, ls="--", c=c)
#     bx.axvline(modnum[i], 0,1,lw=1, ls="--", c=c)
#     ax.plot(modnum, nz, **plot_kwargs)
#     bx.plot(modnum, log_dt, **plot_kwargs)


# def plot_J_tot(hfile, ax="", **plot_kwargs):
#     src, col = getSrcCol(hfile)
#     t = src[:, col.index("star_age")]
#     x = np.log10(t[-1]-t)
#     # x = src[:, col.index("model_number")]
#     J = src[:, col.index("log_total_angular_momentum")]
#     ax.plot(x, J, **plot_kwargs)

# def plot_J_Tc(hfile, ax="", **plot_kwargs):
#     src, col = getSrcCol(hfile)
#     logtc = src[:, col.index("log_center_T")]
#     J = src[:, col.index("log_total_angular_momentum")]
#     ax.plot(logtc, J, **plot_kwargs)


# def plot_J_m(pfile, ax="", show_he_core = False,**plot_kwargs):
#     src, col = getSrcCol(pfile)
#     m = src[:, col.index("mass")]
#     j = src[:, col.index("log_j_rot")]
#     if show_he_core:
#         h1 = src[:, col.index("h1")]
#         k = np.argmin(np.absolute(h1-1e-4))
#         ax.axvline(m[k], 0, 1, ls='--', lw=1, c=c)
#     ax.plot(m, j, **plot_kwargs)


# def plot_Jinside_m(pfile, ax="", show_he_core = False, **plot_kwargs):
#     src, col = getSrcCol(pfile)
#     m = src[:, col.index("mass")]
#     try:
#         j = src[:, col.index("log_J_inside")]
#     except:
#         # column missing, integrate by hand
#         jrot = src[:, col.index("log_j_rot")]
#         j = np.cumsum(jrot)
#     # plot
#     ax.plot(m, j,  c=c, label=label, ls=ls, lw=lw)
#     if show_he_core:
#         h1 = src[:, col.index("h1")]
#         k = np.argmin(np.absolute(h1-1e-4))
#         ax.axvline(m[k], 0, 1, **plot_kwargs)

# def get_mass_rot_res(f):
#     """ given a folder, returns mass and resolution"""
#     model_name = f.split('/')[-2]
#     mass = float(model_name.split('_')[0])
#     initial_rotation = float(model_name.split('rot')[-1].split('_')[0])
#     try:
#         resolution = "$"+model_name.split('_')[2]+model_name.split('_')[3]+"$"
#     except:
#         resolution = "standard"
#     return mass, initial_rotation, resolution


# def get_three_colormaps(z1, z2, z3, Ncolor=10):
#     import matplotlib as mpl
#     """
#     returns three smoothly transitioning colormaps
#     covering the intervals [0, z1], [z1, z2], and [z2, z3]

#     Parameters:
#     ----------
#     z1, z2, z3 : `float` , define the edges of three color regions
#     Ncolor : `int`, optional, number of colors in each segment
#     Returns:
#     --------
#     cmap1, cmap2, cmap3: `mpl.cm` colormaps with corresponding
#     vmin1, vmin2, vmin3, vmax1, vmax2, vmax3: corresponding
#     norm1, norm2, norm3: corresponding norm
#     """
#     # colors that work
#     c1 = 'r'
#     c2 = 'gold'
#     c3 = 'c'
#     c4 = 'fuchsia'

#     cmap1 =  mpl.colors.LinearSegmentedColormap.from_list('cmap1', [c1,c2])
#     vmin1 = 0
#     vmax1 = z1
#     Ncolors1 = Ncolor
#     bounds1 = np.linspace(vmin1, vmax1, Ncolors1)
#     norm1 = mpl.colors.BoundaryNorm(bounds1, cmap1.N)

#     cmap2 = mpl.colors.LinearSegmentedColormap.from_list('cmap2', [c2, c3])
#     vmin2 = vmax1
#     vmax2 = z2
#     Ncolors2 = Ncolor
#     bounds2 = np.linspace(vmin2, vmax2, Ncolors2)
#     norm2 = mpl.colors.BoundaryNorm(bounds2, cmap2.N)

#     cmap3 = mpl.colors.LinearSegmentedColormap.from_list('cmap3', [c3, c4])
#     vmin3 = vmax2
#     vmax3 = z3
#     Ncolors3 = Ncolor
#     bounds3 = np.linspace(vmin3, vmax3, Ncolors3)
#     norm3 = mpl.colors.BoundaryNorm(bounds3, cmap3.N)

#     return(cmap1, vmin1, vmax1, norm1, cmap2, vmin2, vmax2, norm2, cmap3, vmin3, vmax3, norm3)


def get_BE_from_pfile(pfile, alpha_th=1.0, alpha_rot=0.0):
    """Calculates the binding energy profile of the star. See Eq. 6 in
    Dewi & Tauris 2000 (but change sign, binding energy>0 if layer is bound).

    The binding energy is calculated integrating the potential plus a
    fraction `alpha_th` of the internal energy.

    Parameters
    ----------
    pfile    : `MESA profile*.data` file
               assumed to contain the columns mass, energy, radius, and dm
    alpha_th : `float` optional, default 1
               fraction of thermal energy to include in the BE
    Returns
    -------
    BE       : `np.array` binding energy in cgs units

    """
    # from time import time
    # get the data in cgs units
    src, col = getSrcCol(pfile)
    m = src[:, col.index("mass")] * Msun  # g
    dm = src[:, col.index("dm")]  # g
    r = src[:, col.index("radius")] * Rsun_cm  # cm
    u = src[:, col.index("energy")]  # erg
    # calculate local gravitational potential
    psi = -1.0 * G_cgs * np.divide(m, r)  # erg
    if alpha_rot != 0:
        # calculate rotationa energy of spherical shell of
        # mass dm, outer radius r,
        # thickness dr, and rotation frequency omega
        from math import pi
        omega = src[:, col.index("omega")]
        dr = src[:, col.index("dr")]
        I = (2.0/5.0)*dm*((r-dr)**3.0-r**3)
        erot = 0.5*I*omega*omega
    else:
        erot = np.zeros(len(psi))
    # change sign: BE is the energy (>0) to provide to unbind the star
    BE = -1.0 * np.cumsum(np.multiply(psi + alpha_th * u + alpha_rot*erot, dm))
    return np.asarray(BE, dtype=float)


def plot_BE_r(
    pfile,
    ax,
    alpha_th=1.0,
    alpha_rot=0.0,
    top_axis=False,
    mark_he_core=False,
    **plot_kwargs,
):
    """
    plot the binding energy profile as a function of log(r/cm)
    """
    BE = get_BE_from_pfile(pfile, alpha_th=alpha_th, alpha_rot=alpha_rot)
    src, col = getSrcCol(pfile)
    logr = np.log10(src[:, col.index("radius")] * Rsun_cm)
    ax.plot(logr, BE, **plot_kwargs)
    if top_axis:
        """add top axis for mass"""
        ax.tick_params(axis="x", which="both", top=False)
        tx = ax.twiny()
        m = src[:, col.index("mass")]
        xmin, xmax = ax.get_xlim()
        imin = np.argmin(np.absolute(logr - xmin))
        imax = np.argmin(np.absolute(logr - xmax))
        tmin = m[imin]
        tmax = m[imax]
        tx.set_xlim(tmin, tmax)
        # tx.set_xticks([])
        # tx.set_xticklabels([])
        tx.set_xlabel(r"m $[M_\odot]$")
    if mark_he_core:
        m_he_core_01, m_he_core_1, m_he_core_2 = get_He_core_mass_from_pfile(pfile)
        m = src[:, col.index("mass")]
        i_01 = np.argmin(np.absolute(m-m_he_core_01))
        i_1 = np.argmin(np.absolute(m-m_he_core_1))
        i_2 = np.argmin(np.absolute(m-m_he_core_2))
        try:
            color = plot_kwargs['color']
        except:
            color = plot_kwargs['c']
        ax.axvline(logr[i_01], 0, 1, ls="-", lw=2, c=color, zorder=0)
        ax.axvline(logr[i_1], 0, 1, ls="--", lw=2, c=color, zorder=0)
        ax.axvline(logr[i_2], 0, 1, ls="-.", lw=2, c=color, zorder=0)


def plot_BE_m(pfile, ax, alpha_th=1.0, **plot_kwargs):
    """
    plot the binding energy profile as a function of mass
    """
    BE = get_BE_from_pfile(pfile, alpha_th=alpha_th)
    src, col = getSrcCol(pfile)
    # logr = np.log10(src[:, col.index("radius")] * Rsun_cm)
    m = src[:, col.index("mass")]
    ax.plot(m, BE, **plot_kwargs)


def get_ratio_BE(pfile1, pfile2, alpha_th, alpha_rot):
    """Ratio of pfile2/pfile1 on the x-coordinates in pfile1.

    Ratio of the binding energies (0<alpha_th<=1
    sets the fractional contribution of the internal energy -- 0<alpha_rot<=1
    sets the fractional contribution of the rotational energy) or
    gravitational energy (alpha_th=0, alpha_rot=0).

    Uses scipy.interpolate.interp1d to find a function that can then
    be evaluated on the other grid. The interpolation is done in
    fractional mass coordinate q=m/M.

    Parameters:
    -----------
    pfile1   :   `string` absolute path of MESA profile
    pfile2   :   `string` absolute path of MESA profile
    alpha_th :   `float`, optional fraction of thermal energy that is
                 considered in the binding energy. By default include the internal energy.
    alpha_rot:   `float`, optional fraction of rotational energy that is
                 considered in the binding energy. By default do not the rotational energy.

    Returns:
    -------
    ratio2_to_1: `np.array`, array of BE from pfile2/BE from pfile1 (or np.nan for the domain
                  of pfile2 outside the domain of pfile1 in mass)

    """
    # get data pfile1
    src1, col1 = getSrcCol(pfile1)
    BE1 = get_BE_from_pfile(pfile1, alpha_th=alpha_th, alpha_rot=alpha_rot)
    r1 = src1[:, col1.index("radius")]
    m1 = src1[:, col1.index("mass")]
    q1 = m1 / max(m1)
    # get data pfile2
    src2, col2 = getSrcCol(pfile2)
    BE2 = get_BE_from_pfile(pfile2, alpha_th=alpha_th, alpha_rot=alpha_rot)
    r2 = src2[:, col2.index("radius")]
    m2 = src2[:, col2.index("mass")]
    q2 = m2 / max(m2)
    # interpolate pfile2 on pfile1
    # MESA arrays are from the surface inwards, to have a monotonically increasing
    # mass coordinate, flip the arrays to interpolate
    x1 = q1[::-1]
    x2 = q2[::-1]
    y1 = BE1[::-1]
    y2 = BE2[::-1]
    interp2 = interp1d(
        x2, y2, assume_sorted=True, bounds_error=False, fill_value="extrapolate"
    )
    y2_interp = interp2(x1)
    ratio2_to_1 = y2_interp / y1
    # flip back to MESA convention (index 0 == surface)
    ratio2_to_1 = ratio2_to_1[::-1]
    return ratio2_to_1


def plot_ratio_BE_r(pfile1, pfile2, ax, alpha_th=1.0, alpha_rot=0.0,  **plot_kwargs):
    """
    Ratio of pfile2/pfile1 on the x-coordinates in pfile1.

    Ratio of the binding energies (0<alpha_th<=1
    sets the fractional contribution of the internal energy -- 0<alpha_rot<=1
    sets the fractional contribution of the rotational energy) or
    gravitational energy (alpha_th=0, alpha_rot=0).

    Uses scipy.interpolate.interp1d to find a function that can then
    be evaluated on the other grid. The interpolation is done in mass coordinate.
    Extrapolated values are set to nan (very center and very outermost layers typically).

    Parameters:
    -----------
    pfile1   :   `string` absolute path of MESA profile
    pfile2   :   `string` absolute path of MESA profile
    ax       :   `matplotlib.axes._subplots.AxesSubplot`
    alpha_th :   `float`, optional fraction of thermal energy that is
                considered in the binding energy. By default include the internal energy.
    alpha_rot:   `float`, optional fraction of rotational energy that is
                 considered in the binding energy. By default do not the rotational energy.


    Returns:
    -------
    ratio2_to_1: `np.array`, array of BE from pfile2/BE from pfile1 (or np.nan for the domain
                  of pfile2 outside the domain of pfile1 in mass)
    """
    src1, col1 = getSrcCol(pfile1)
    ratio2_to_1 = get_ratio_BE(pfile1, pfile2, alpha_th, alpha_rot)
    r1 = src1[:, col1.index("radius")]
    ax.plot(np.log10(r1 * Rsun_cm), ratio2_to_1, **plot_kwargs)
    return ratio2_to_1


def plot_ratio_BE_m(pfile1, pfile2, ax, alpha_th=1, **plot_kwargs):
    """
    Ratio of pfile2/pfile1 on the x-coordinates in pfile1.

    plots the ratio of the binding energies (0<alpha_th<=1
    sets the fractional contribution of the internal energy) or
    gravitational energy (alpha_th=0).

    Uses scipy.interpolate.interp1d to find a function that can then
    be evaluated on the other grid. The interpolation is done in mass coordinate.
    Extrapolated values are set to nan (very center and very outermost layers typically).

    Parameters:
    -----------
    pfile1  :   `string` absolute path of MESA profile
    pfile2  :   `string` absolute path of MESA profile
    ax      :   `matplotlib.axes._subplots.AxesSubplot`
    alpha_th:   `float`, optional fraction of thermal energy that is
                considered in the binding energy. By default include the internal energy.

    Returns:
    -------
    ratio2_to_1: `np.array`, array of BE from pfile2/BE from pfile1 (or np.nan for the domain
                  of pfile2 outside the domain of pfile1 in mass)
    """
    # get data pfile1
    src1, col1 = getSrcCol(pfile1)
    BE1 = get_BE_from_pfile(pfile1, alpha_th=alpha_th)
    m1 = src1[:, col1.index("mass")]
    # get data pfile2
    src2, col2 = getSrcCol(pfile2)
    BE2 = get_BE_from_pfile(pfile2, alpha_th=alpha_th)
    m2 = src2[:, col2.index("mass")]
    # interpolate pfile2 on pfile1
    # do the interpolation in mass coordinate and then plot in radius
    # MESA arrays are from the surface inwards, to have a monotonically increasing
    # mass coordinate, flip the arrays to interpolate
    x1 = m1[::-1]
    x2 = m2[::-1]
    y1 = BE1[::-1]
    y2 = BE2[::-1]
    # restrict the interpolation to the range of overlap
    xmin = max(min(x1), min(x2))
    xmax = min(max(x1), max(x2))
    ind = (x2 >= xmin) & (x2 <= xmax)
    interp2 = interp1d(x2[ind], y2[ind], assume_sorted=True, bounds_error=False)
    y2_interp = interp2(x1)
    ratio2_to_1 = y2_interp / y1
    ax.plot(m1, ratio2_to_1, **plot_kwargs)
    return ratio2_to_1


# def plot_ratio_Eint_r(pfile1, pfile2, ax, **plot_kwargs):
#     """plots the ratio of the internal energies. See also plot_ratio_BE_r.

#     Parameters:
#     -----------
#     pfile1 : `string` absolute path of MESA profile
#     pfile2 : `string` absolute path of MESA profile
#     ax     : `matplotlib.axes._subplots.AxesSubplot`

#     Returns:
#     -------
#     ratio2_to_1: `np.array`, array of BE from pfile2/BE from pfile1 (or np.nan for the domain
#                   of pfile2 outside the domain of pfile1 in mass)

#     """
#     # get data pfile1
#     src1, col1 = getSrcCol(pfile1)
#     e1 = src1[:, col1.index("energy")]
#     logr1 = np.log10(src1[:, col1.index("radius")] * Rsun_cm)
#     m1 = src[:, col.index("mass")]
#     # get data pfile2
#     src2, col2 = getSrcCol(pfile2)
#     e2 = src2[:, col2.index("energy")]
#     logr2 = np.log10(src2[:, col2.index("radius")] * Rsun_cm)
#     m2 = src[:, col.index("mass")]
#     x1 = m1
#     x2 = m2
#     y1 = e1
#     y2 = e2
#     #
#     # MESA arrays are from the surface inwards, to have a monotonically increasing
#     # mass coordinate, flip the arrays to interpolate
#     interp2 = interp1d(x2[::-1], y2[::-1], assume_sorted=True, bounds_error=False)
#     y2_interp = interp2(x1)
#     ratio2_to_1 = y2_interp / y1
#     # print("2_to_1", np.nanmax(ratio2_to_1), np.nanmin(ratio2_to_1))
#     ax.plot(np.log10(r1 * Rsun_cm), ratio2_to_1, **plot_kwargs)
#     return ratio2_to_1


# def plot_ratio_Eint_m(pfile1, pfile2, ax, **plot_kwargs):
#     """plots the ratio of the internal energies. See also plot_ratio_BE_r.

#     Parameters:
#     -----------
#     pfile1 : `string` absolute path of MESA profile
#     pfile2 : `string` absolute path of MESA profile
#     ax     : `matplotlib.axes._subplots.AxesSubplot`

#     Returns:
#     -------
#     ratio2_to_1: `np.array`, array of BE from pfile2/BE from pfile1 (or np.nan for the domain
#                   of pfile2 outside the domain of pfile1 in mass)

#     """
#     # get data pfile1
#     src1, col1 = getSrcCol(pfile1)
#     e1 = src1[:, col1.index("energy")]
#     m1 = src[:, col.index("mass")]
#     # get data pfile2
#     src2, col2 = getSrcCol(pfile2)
#     e2 = src2[:, col2.index("energy")]
#     m2 = src[:, col.index("mass")]
#     x1 = m1
#     x2 = m2
#     y1 = e1
#     y2 = e2
#     #
#     # MESA arrays are from the surface inwards, to have a monotonically increasing
#     # mass coordinate, flip the arrays to interpolate
#     interp2 = interp1d(x2[::-1], y2[::-1], assume_sorted=True, bounds_error=False)
#     y2_interp = interp2(x1)
#     ratio2_to_1 = y2_interp / y1
#     # print("2_to_1", np.nanmax(ratio2_to_1), np.nanmin(ratio2_to_1))
#     ax.plot(m1, ratio2_to_1, **plot_kwargs)
#     return ratio2_to_1


def get_ax_from_pfile(pfile, axes_list):
    """
    Parameters:
    ----------
    pfile : `string` absolute path to MESA profile file, must contain one of the strings below
    axes_list : `list` of mpl.axes objects
    Returns:
    -------
    one of the elements of the list
    """
    if "100Rsun" in pfile:
        return axes_list[0]
    elif "200Rsun" in pfile:
        return axes_list[1]
    elif "300Rsun" in pfile:
        return axes_list[2]
    elif "500Rsun" in pfile:
        return axes_list[3]
    elif "1000Rsun" in pfile:
        return axes_list[4]
    else:
        print(colored(pfile, "red"))
        raise ValueError(pfile + " not a known radius checkpoint")


# def get_radius(pfile):
#     """ returns outer radius of profile file in Rsun units """
#     src, col = getSrcCol(pfile)
#     radius = src[:, col.index("radius")]
#     # print(max(radius))
#     return max(radius)


# def radius_sort(pfile):
#     """ returns radius from profile name"""
#     return int(pfile.split('/')[-1].split('.data')[0].rstrip("Rsun"))


def get_lambda_profile(pfile, alpha_th=1):
    """Calculates lambda_CE according to Eq. 1 of Dewi & Tauris 2000 or Eq. 9 of de Kool 1990
    (N.B.: Ivanova et al. 2013, fig. 5 y-axis label is 1/lambda not lambda)

    Parameters:
    --------
    pfile    : `string` path MESA profile*.data file
               assumed to contain the columns mass, energy, radius, and dm
    alpha_th : `float` optional
               fraction of thermal energy to include in the BE
    Returns
    -------
    L        : `np.array`, lambda values for each mesh point (same grid as input pfile)
    """
    # get the data in cgs units, load from center outwards
    src, col = getSrcCol(pfile)
    m = src[:, col.index("mass")]  # Msun
    M = max(m)
    # the radius is binary_separation*roche_lobe_radius, so fixed to the surface
    R = src[0, col.index("radius")]  # Rsun
    BE = get_BE_from_pfile(pfile, alpha_th=alpha_th)
    L = (G_cgs * M * (M - m) / R) / BE
    # add units back
    L *= Msun * Msun / Rsun_cm
    # print(f"{BE[i]:.5e}")
    # print(m[i]/Msun, (M-m[-1])/Msun, BE[-1])
    # print(L[i])
    # print(G_cgs)
    # print("-----------")
    return np.asarray(L, dtype=float)


def interpolate_BE_env(x1, x2, y2, offset=1.0):
    """
    Interpolates y2(x2) to get y2(x1)

    Parameters:
    ----------
    x1: `np.array`, independent coords we want
    y2, x2: `np.array`, need the same length
    offset : `float`, offset to make smaller numbers for interpolation

    Returns:
    --------
    y2_interp == y2(x1)
    """
    interpolator = interp1d(
        x2, y2 / offset, assume_sorted=True, bounds_error=False, fill_value=np.nan
    )
    y2_interp = interpolator(x1) * offset
    return np.array(y2_interp, dtype=float)


def plot_lambda_mass(pfile, ax, alpha_th=1, **plot_kwargs):
    l = get_lambda_profile(pfile, alpha_th=alpha_th)
    src, col = getSrcCol(pfile)
    m = src[:, col.index("mass")]
    ax.plot(m, l, **plot_kwargs)


def plot_lambda_r(pfile, ax, alpha_th=1, **plot_kwargs):
    l = get_lambda_profile(pfile, alpha_th=alpha_th)
    src, col = getSrcCol(pfile)
    r = src[:, col.index("radius")]
    ax.plot(r, l, **plot_kwargs)


def get_He_core_mass_from_pfile(pfile):
    """
    Returns the He core mass based on the default MESA
    definition
    """
    src, col = getSrcCol(pfile)
    x = src[:, col.index("h1")]
    y = src[:, col.index("he4")]
    m = src[:, col.index("mass")]
    Hecore = (y >= 0.1) & (x <= 0.01)
    m_he_core_01 = max(m[Hecore])
    Hecore = (y >= 0.1) & (x <= 0.1)
    m_he_core_1 = max(m[Hecore])
    Hecore = (y >= 0.1) & (x <= 0.2)
    m_he_core_2 = max(m[Hecore])
    return m_he_core_01, m_he_core_1, m_he_core_2


def plot_lambda_at_one_radius(
    ax, string, grid_folders, accretor=None, nonrot=None, plot_func=plot_lambda_mass
):
    """make plots of the lambda profiles comparing models

    Parameters:
    ----------
    ax :          `mpl.axes` where to plot
    string:       `str` profile name to be looked into the folders LOGS/LOGS2
    grid_folders: `str` path to grid of engineered models workdirectory
    accretor:     `str` or None, optional, path to the LOGS2 folder of the accretor
    nonrot:       `str` or None, optional, path to the work directory of a single star normal model
    plot_func  :  `python function`, can be `plot_lambda_mass` or
                  `plot_lambda_r` depending on if you want to plot
                  as a function of mass coordinate or radius
    legend:       `bool` whether to show the legend or not
    """
    colors = plt.cm.viridis(np.linspace(0, 1, len(grid_folders)))

    # engineered stars
    for f in grid_folders:
        c = colors[grid_folders.index(f)]
        pfile = glob.glob(f + string)[0]
        plot_func(pfile, ax, alpha_th=1, c=c, lw=2, label=f.split("/")[-2])

    if accretor:
        pfile = glob.glob(accretor + "/" + string)[0]
        plot_func(pfile, ax, alpha_th=1, c="orange", label="accretor")
        m_he_core_01, m_he_core_1, m_he_core_2 = get_He_core_mass_from_pfile(pfile)
        ax.axvline(m_he_core_01, 0, 1, ls="-", lw=2, c="orange", zorder=0)
        ax.axvline(m_he_core_1, 0, 1, ls="--", lw=2, c="orange", zorder=0)
        ax.axvline(m_he_core_2, 0, 1, ls="-.", lw=2, c="orange", zorder=0)

    if nonrot:
        pfile = glob.glob(nonrot + string)[0]
        plot_func(pfile, ax, alpha_th=1, c="r", label="single")
        m_he_core_01, m_he_core_1, m_he_core_2 = get_He_core_mass_from_pfile(pfile)
        ax.axvline(m_he_core_01, 0, 1, ls="-", lw=2, c="red", zorder=0)
        ax.axvline(m_he_core_1, 0, 1, ls="--", lw=2, c="red", zorder=0)
        ax.axvline(m_he_core_2, 0, 1, ls="-.", lw=2, c="red", zorder=0)
