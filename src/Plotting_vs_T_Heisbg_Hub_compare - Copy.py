import h5py
import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from mpl_toolkits import mplot3d
from itertools import combinations
from itertools import permutations

# file name to load
file_name = (
    "./results/M_vs_T/Heisbg_few_T_sweep_h1_3_result_08_Dec_2022_12_39_nwynn=3.jld"
)

# which slice of h to plot
h_inds = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10]

# parameters for plotting
# This is to optionally skip very low temperatures data in plotting
# as they can get wild. Tlimit sets the index of the lowest temperature ploted.
Tlimit = 0
hlimit = 0
Tupperlimit = -1
hupperlimit = 0


# load data
# cd to the folder containing the script
abspath = os.path.abspath(__file__)
dname = os.path.dirname(abspath)
os.chdir(dname)

# offload parameters and result from data
data = h5py.File(file_name, "r")
M_order8 = np.array(data["M11"].value)
M_order9 = np.array(data["M12"].value)
S_order8 = np.array(data["S11"].value)
S_order9 = np.array(data["S12"].value)
M_Euler_order8 = np.array(data["M_Euler11"].value)
M_Euler_order9 = np.array(data["M_Euler12"].value)
S_Euler_order8 = np.array(data["S_Euler11"].value)
S_Euler_order9 = np.array(data["S_Euler12"].value)
M_raw_order8 = np.array(data["M_raw11"].value)
M_raw_order9 = np.array(data["M_raw12"].value)
S_raw_order8 = np.array(data["S_raw11"].value)
S_raw_order9 = np.array(data["S_raw12"].value)
Temps = np.array(data["Temps"].value)
hs = np.array(data["hs"].value)


def Mplotting():

    # Number of temperatures
    NT = len(Temps)
    # Number of magnetic field h
    Nh = len(hs)

    # wynn sum data

    tot_h_nums = len(h_inds)

    # initialize the plot
    plt.figure(figsize=(10, 8))
    plt.rc("axes", labelsize=15)
    plt.rc("xtick", labelsize=15)
    plt.rc("ytick", labelsize=15)
    plt.rc("legend", fontsize=6)

    colors = ["g", "r", "b", "c", "k", "m", "y"]

    # plt.rcParams['font.size']= 15
    # Energy
    ax = plt.subplot(111)

    # Magnetization NLCE
    plt.subplot(111)
    for h_num, h_index in enumerate(h_inds):
        plt.plot(
            Temps[Tlimit:Tupperlimit],
            M_order8[h_index, Tlimit:Tupperlimit],
            "-",
            color=cm.rainbow(h_num / tot_h_nums),
            label="Wynn sum, order11" + ", h=" + str(hs[h_index]),
        )
        plt.plot(
            Temps[Tlimit:Tupperlimit],
            M_order9[h_index, Tlimit:Tupperlimit],
            ".",
            color=cm.rainbow(h_num / tot_h_nums),
            label="Wynn sum, order12" + ", h=" + str(hs[h_index]),
        )
    plt.axhline(1 / 6)
    plt.axhline(1 / 4)
    plt.axhline(1 / 3)
    ax.set_xlabel("T")
    ax.set_ylabel("M")
    plt.ylim(0.0, 0.52)
    plt.legend(loc="upper right", frameon=False)

    # ======================== 3d plot of S ========================
    """
    plt.figure(figsize=(10, 8))
    plt.rc("axes", labelsize=15)
    plt.rc("xtick", labelsize=15)
    plt.rc("ytick", labelsize=15)
    plt.rc("legend", fontsize=15)

    ax = plt.axes(projection="3d")
    h, T = np.meshgrid(hs[hlimit:], Temps[Tlimit:])
    ax.contour3D(h, T, S[Nmax - 1, :, :], 200, cmap="binary")
    ax.set_xlabel("h")
    ax.set_ylabel("T")
    ax.set_zlabel("S")

    """
    """
    plt.contourf(
        Main.hs[hlimit:], Temp[Tlimit:], entropy, 500, label="%s" % N, cmap="RdGy",
    )
    """
    """

    # plt.colorbar()
    plt.ylabel("T")
    plt.xlabel("h")
    # plt.xlim(0.1,10)
    # plt.ylim(0.0, 1.5)
    # plt.xscale("log")
    # plt.legend(loc='lower right')
    # plt.tight_layout()
    # plt.subplots_adjust(right = 1.1)
    # plt.show()


    plt.savefig(
        "S_3d_Hubbard_triangular_U=" + str(U) + "_order_" + str(Nmax) + ".png",
        format="png",
        dpi=600,
        transparent=False,
    )"""


def Splotting():
    # Number of temperatures
    NT = len(Temps)
    # Number of magnetic field h
    Nh = len(hs)

    # initialize the plot
    plt.figure(figsize=(10, 8))
    plt.rc("axes", labelsize=15)
    plt.rc("xtick", labelsize=15)
    plt.rc("ytick", labelsize=15)
    plt.rc("legend", fontsize=6)

    colors = ["g", "r", "b", "c", "k", "m", "y"]

    # plt.rcParams['font.size']= 15
    # Energy
    ax = plt.subplot(111)

    # Magnetization NLCE
    plt.subplot(111)

    # wynn sum data
    tot_h_nums = len(h_inds)
    for h_num, h_index in enumerate(h_inds):
        plt.plot(
            Temps[Tlimit:Tupperlimit],
            S_order8[h_index, Tlimit:Tupperlimit],
            "-",
            color=cm.rainbow(h_num / tot_h_nums),
            label="wynn sum, order11" + ", h=" + str(hs[h_index]),
        )
        plt.plot(
            Temps[Tlimit:Tupperlimit],
            S_order9[h_index, Tlimit:Tupperlimit],
            ".",
            color=cm.rainbow(h_num / tot_h_nums),
            label="wynn sum, order12" + ", h=" + str(hs[h_index]),
        )
    ax.set_xlabel("T")
    ax.set_ylabel("S")
    plt.ylim(-0.05, 0.7)
    plt.legend(loc="upper right", frameon=False)


Mplotting()

plt.savefig(
    "./plots/M_vs_T_medium_h" +
    # "h=" + str(hs[h_index]) + "_order_" +
    "_Wynn_11&12_3cycles" + ".png",
    format="png",
    dpi=600,
    transparent=False,
)

Splotting()

plt.savefig(
    "./plots/S_vs_T_medium_h" +
    # "h=" + str(hs[h_index]) + "_order_" +
    "_Wynn_11&12_3cycles" + ".png",
    format="png",
    dpi=600,
    transparent=False,
)
