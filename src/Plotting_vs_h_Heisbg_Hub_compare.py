import h5py
import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from mpl_toolkits import mplot3d
from itertools import combinations
from itertools import permutations

############################
##  Do Resummation here   ##
############################
# don't use this part yet! It's still bugged!!!!!
new_resum = False
# only run this part if do new resummation.
if new_resum:
    import julia
    from julia import Main

    Main.include("./resummation.jl")
    Main.Nmax = 12
    Main.raw_sum_order = 11
    Main.nwynn = 4

    # run the resummation to the desired order
    Main.eval(
        "out_name = print_data(Nmax=Nmax, raw_sum_order=raw_sum_order, nwynn=nwynn)"
    )

    # file name to load
    file_name = Main.out_name
# "./results/Heisbg_few_T_sweep_h_result_30_Nov_2022_23_58.jld"

###############################
##  Do plotting below        ##
###############################

# Only run this part if don't do new resummation.
# Load existing result file in that case.
if not new_resum:
    file_name = (
        "./results/Heisbg_few_T_sweep_h_result_08_Dec_2022_15_22_nwynn5_Euler4.jld"
    )

# which slice of h to plot
T_inds = [2, 3, 4, 5, 7]
# [0, 1, 2, 3, 4, 5, 6, 7, 8]

# parameters for plotting
# This is to optionally skip very low temperatures data in plotting
# as they can get wild. hlimit sets the index of the lowest temperature ploted.
hlimit = 0
hlimit = 0
hupperlimit = 0
hupperlimit = -1


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


def M_vs_h_plotting():

    # Number of temperatures
    NT = len(Temps)
    # Number of magnetic field h
    Nh = len(hs)

    # wynn sum data

    tot_T_nums = len(T_inds)

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
    for T_num, T_index in enumerate(T_inds):
        print(
            "Chi at T="
            + str(Temps[T_index])
            + ", h=0 is "
            + str(M_order9[1, T_index] / hs[1])
        )
        plt.plot(
            hs[hlimit:hupperlimit],
            M_order8[hlimit:hupperlimit, T_index],
            "--",
            color=cm.rainbow(T_num / tot_T_nums),
            label="wynn sum, order11" + ", T=" + str(Temps[T_index]),
        )
        plt.plot(
            hs[hlimit:hupperlimit],
            M_order9[hlimit:hupperlimit, T_index],
            "-",
            color=cm.rainbow(T_num / tot_T_nums),
            label="wynn sum, order12" + ", T=" + str(Temps[T_index]),
        )
        """
        plt.plot(
            hs[hlimit:hupperlimit],
            M_raw_order9[hlimit:hupperlimit, T_index],
            "-.",
            color=cm.rainbow(T_num / tot_T_nums),
            label="raw sum, order12" + ", T=" + str(Temps[T_index]),
        )
        plt.plot(
            hs[hlimit:hupperlimit],
            M_raw_order8[hlimit:hupperlimit, T_index],
            "--",
            color=cm.rainbow(T_num / tot_T_nums),
            label="raw sum, order11" + ", T=" + str(Temps[T_index]),
        )
        """
        plt.plot(
            hs[hlimit:hupperlimit],
            M_Euler_order9[hlimit:hupperlimit, T_index],
            ".",
            color=cm.rainbow(T_num / tot_T_nums),
            label="Euler sum, order12" + ", T=" + str(Temps[T_index]),
        )
        plt.axhline(1 / 6)
        plt.axhline(1 / 4)
        plt.axhline(1 / 3)
    ax.set_xlabel("h")
    ax.set_ylabel("M")
    plt.ylim(0.0, 0.52)
    plt.legend(loc="lower right", frameon=False)

    # ======================== 3d plot of S ========================
    """
    plt.figure(figsize=(10, 8))
    plt.rc("axes", labelsize=15)
    plt.rc("xtick", labelsize=15)
    plt.rc("ytick", labelsize=15)
    plt.rc("legend", fontsize=15)

    ax = plt.axes(projection="3d")
    h, T = np.meshgrid(hs[hlimit:], Temps[hlimit:])
    ax.contour3D(h, T, S[Nmax - 1, :, :], 200, cmap="binary")
    ax.set_xlabel("h")
    ax.set_ylabel("T")
    ax.set_zlabel("S")

    """
    """
    plt.contourf(
        Main.hs[hlimit:], Temp[hlimit:], entropy, 500, label="%s" % N, cmap="RdGy",
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


def S_vs_h_plotting():
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
    tot_T_nums = len(T_inds)
    for T_num, T_index in enumerate(T_inds):
        plt.plot(
            hs[hlimit:hupperlimit],
            S_order8[hlimit:hupperlimit, T_index],
            "--",
            color=cm.rainbow(T_num / tot_T_nums),
            label="wynn sum, order11" + ", T=" + str(Temps[T_index]),
        )
        plt.plot(
            hs[hlimit:hupperlimit],
            S_order9[hlimit:hupperlimit, T_index],
            "-",
            color=cm.rainbow(T_num / tot_T_nums),
            label="wynn sum, order12" + ", T=" + str(Temps[T_index]),
        )
        """
        plt.plot(
            hs[hlimit:hupperlimit],
            S_raw_order9[hlimit:hupperlimit, T_index],
            "-",
            color=cm.rainbow(T_num / tot_T_nums),
            label="raw sum, order9" + ", T=" + str(Temps[T_index]),
        )
        plt.plot(
            hs[hlimit:hupperlimit],
            S_raw_order8[hlimit:hupperlimit, T_index],
            "-.",
            color=cm.rainbow(T_num / tot_T_nums),
            label="raw sum, order8" + ", T=" + str(Temps[T_index]),
        )
        plt.plot(
            hs[hlimit:hupperlimit],
            S_Euler_order9[hlimit:hupperlimit, T_index],
            ".",
            color=cm.rainbow(T_num / tot_T_nums),
            label="Euler sum, order9" + ", T=" + str(Temps[T_index]),
        )"""
    ax.set_xlabel("h")
    ax.set_ylabel("S")
    plt.ylim(-0.05, 0.7)
    plt.legend(loc="upper right", frameon=False)


M_vs_h_plotting()

plt.savefig(
    "./plots/M_vs_h_Heisbg_triangular_multi_T" +
    # "h=" + str(hs[T_index]) + "_order_" +
    "_Wynn&raw&Euler_8&9" + ".png",
    format="png",
    dpi=600,
    transparent=False,
)


S_vs_h_plotting()

plt.savefig(
    "./plots/S_vs_h_Heisbg_triangular_multi_T" +
    # "h=" + str(hs[T_index]) + "_order_" +
    "_Wynn&raw&Euler_8&9" + ".png",
    format="png",
    dpi=600,
    transparent=False,
)

