import h5py
import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from mpl_toolkits import mplot3d
from itertools import combinations
from itertools import permutations

# file name to load
file_name = "./results/Heisbg_exp_compare/Heisbg_exp_compare_20_Dec_2022_20_39.jld"

exp_data = np.loadtxt(
    "../../expriment_fitting/exp_data/M_T.dat", skiprows=2, max_rows=2578
)

# which slice of h to plot
h_inds = [0, 1, 2, 3, 4]

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
M_order9 = np.array(data["M12"].value)
S_order9 = np.array(data["S12"].value)
M_Euler_order9 = np.array(data["M_Euler12"].value)
S_Euler_order9 = np.array(data["S_Euler12"].value)
M_raw_order9 = np.array(data["M_raw12"].value)
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
            exp_data[:, 2 * h_index],
            exp_data[:, 2 * h_index + 1],
            "-",
            color=cm.rainbow(h_num / tot_h_nums),
            label="exp" + ", h=" + str(hs[h_index]),
        )
        plt.plot(
            Temps[Tlimit:Tupperlimit],
            M_order9[h_index, Tlimit:Tupperlimit] * 10 ** 3,
            ".",
            color=cm.rainbow(h_num / tot_h_nums),
            label="Wynn sum, order12" + ", h=" + str(hs[h_index]),
        )
    # plt.axhline(1 / 6)
    # plt.axhline(1 / 4)
    # plt.axhline(1 / 3)
    ax.set_xlabel("T")
    ax.set_ylabel(r"$M$ (emu/mol_Cu)")
    plt.xlim(0.0, 50)
    # plt.ylim(0.0, 0.52)
    plt.legend(loc="upper right", frameon=False)


Mplotting()

plt.savefig(
    "./plots/M_vs_T_exp_compare" +
    # "h=" + str(hs[h_index]) + "_order_" +
    "_Wynn_12_4cycles" + ".png",
    format="png",
    dpi=600,
    transparent=False,
)

