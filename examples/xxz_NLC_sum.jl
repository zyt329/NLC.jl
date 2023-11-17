using JLD2
using NLC
using Printf


"""
    Max_num_clusters defines maximum number of clusters in the highest order. Setting it high may result in insufficient memory.

    Input: only folders that contain simulation info. Parameters can be read from the folders.
"""
function run_NLC_sum(; Nmax)

    #================================#
    #= DEFINE SIMULATION PARAMETERS =#
    #================================#

    # don't need to define since it would be automatically passed in 

    # ==============================
    # ====  make/find folders   ====
    # ==============================

    # path to diagonalization folder, in case thermal data was damaged
    diag_folder_path = "./test_diagonization_data_J_z[J_xy1.0000-1"

    # path to cluster info
    clusters_info_path = "../cluster_info/triangle/"

    # pass in simulation folders
    simulation_folder_full_path = "./test_simulation_J_z[J_xy1.0000-1"

    # ==============================
    # ==== thermally average  ======
    # ==============================

    # do the thermal averages and find out multiplicities of clusters
    thermal_avg_folder_full_path = simulation_folder_full_path * "/thermal_avg_data/"

    # ===========================
    # ==== do the NLC sum  ======
    # ===========================


    O = NLC_sum(Nmax=Nmax, clusters_info_path=clusters_info_path, simulation_folder_full_path=simulation_folder_full_path, thermal_avg_folder_full_path=thermal_avg_folder_full_path, diag_folder_path=diag_folder_path)


end

@time run_NLC_sum(Nmax=9)

