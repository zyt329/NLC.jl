using NLC
using Printf

function run_diagonalization(ARGS)

    #================================#
    #= DEFINE SIMULATION PARAMETERS =#
    #================================#

    # xy interaction
    J_xy = 1.0

    # z interaction
    J_z = parse(Float64, ARGS[1])

    # g factor
    g = 1.0

    # ===================================================
    # ==== make folder to hold diagonalization data  ====
    # ===================================================

    diag_folder_prefix = @sprintf "diagonization_data_J_z[J_xy%.4f" (J_z / J_xy)

    diag_folder_full_path = make_indexed_folder(folder_prefix=diag_folder_prefix)

    # ==============================
    # ==== start diagonalizing  ====
    # ==============================

    diagonalize_all_clusters_xxz(; J_xy=J_xy, J_z=J_z, Nmax=12, clusters_info_path="/nfs/home/zyt329/Heisbg_triang/xxz/NLC_cluster_info/", diag_folder_path=diag_folder_full_path)

end


"""
    Max_num_clusters defines maximum number of clusters in the highest order. Setting it high may result in insufficient memory.
"""
function run_NLCE(ARGS; Nmax, Max_num_clusters)

    #================================#
    #= DEFINE SIMULATION PARAMETERS =#
    #================================#

    # temperature values
    Temps = vcat(range(0.0, 5, length=51))
    NT = length(Temps)

    # field values
    hs = vcat(range(0.0, 2.0, length=5))
    Nh = length(hs)

    # Jz/Jxy ratio
    J_z_by_J_xy = parse(Float64, ARGS[1])

    # auxiliary variable θ to keep jz and jxy normalized such that Jz^2 + J_xy^2 = 1
    θ = atan(J_z_by_J_xy)

    # xy interaction
    J_xy = cos(θ)

    # z interaction
    J_z = sin(θ)

    # g factor
    g = 1.0

    # ==============================
    # ====  make/find folders   ====
    # ==============================

    # path to cluster info
    clusters_info_path = "/nfs/home/zyt329/Heisbg_triang/xxz/NLC_cluster_info/"

    # find folder with diagonalization data
    diag_folder_full_path = @sprintf "diagonization_data_J_z[J_xy%.4f-1" (J_z / J_xy)

    # make simulation folders
    simulation_folder_prefix = @sprintf "simulation_J_z[J_xy%.4f" (J_z / J_xy)

    simulation_folder_full_path = make_indexed_folder(folder_prefix=simulation_folder_prefix, folder_path=".")

    # ==============================
    # ==== thermally average  ======
    # ==============================

    # do the thermal averages and find out multiplicities of clusters
    multi, thermal_avg_folder_full_path = thermal_avg_all_clusters_xxz(Nmax=Nmax, Max_num_clusters=Max_num_clusters, J_xy=J_xy, J_z=J_z, g=1, Temps=Temps, hs=hs, diag_folder_full_path=diag_folder_full_path, simulation_folder_full_path=simulation_folder_full_path, clusters_info_path=clusters_info_path)

    # ===========================
    # ==== do the NLC sum  ======
    # ===========================

    O = NLC_sum(Nmax=Nmax, Max_num_clusters=Max_num_clusters, J_xy=J_xy, J_z=J_z, g=1, Temps=Temps, hs=hs, multi=multi, clusters_info_path=clusters_info_path, simulation_folder_full_path=simulation_folder_full_path, thermal_avg_folder_full_path=thermal_avg_folder_full_path)

end

run_NLCE(ARGS, Nmax=9, Max_num_clusters=1500)
