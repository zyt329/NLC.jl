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

    #=
    # temperature values
    Temps = vcat(range(0.0, 1.0, length=51))
    NT = length(Temps)

    # field values
    hs = vcat(range(0.0, 6.0, length=41))
    Nh = length(hs)

    """# Jz/Jxy ratio
    J_z_by_J_xy = J_z_by_J_xy

    # auxiliary variable θ to keep jz and jxy normalized such that Jz^2 + J_xy^2 = 1
    θ = atan(J_z_by_J_xy)

    # xy interaction
    J_xy = cos(θ)

    # z interaction
    J_z = sin(θ)"""

    # g factor
    g = 1.0
    =#

    # ==============================
    # ====  make/find folders   ====
    # ==============================


    # path to cluster info
    clusters_info_path = "/gpfs/scratch/yzhang/Research/xxz/NLC_clusters_info_JSON/triangle_clusters_no_2_order/"

    # make simulation folders
    #simulation_folder_prefix = @sprintf "hashtag_simulation_J_z[J_xy%.4f" (J_z / J_xy)

    simulation_folder_full_path = "/gpfs/scratch/yzhang/Research/xxz/runs/triangle_based_no2_Ising_simulation-1/"

    #make_indexed_folder(folder_prefix=simulation_folder_prefix, folder_path="/nfs/home/zyt329/Research/xxz/runs/")

    # ==============================
    # ==== thermally average  ======
    # ==============================

    # do the thermal averages and find out multiplicities of clusters
    thermal_avg_folder_full_path = simulation_folder_full_path * "/thermal_avg_data/"

    # ===========================
    # ==== do the NLC sum  ======
    # ===========================


    O = NLC_sum_Ising(Nmax=Nmax, clusters_info_path=clusters_info_path, simulation_folder_full_path=simulation_folder_full_path, thermal_avg_folder_full_path=thermal_avg_folder_full_path)


end



for J_z_by_J_xy in [1]
    # run_diagonalization(J_z=Float64(J_z))
    @time run_NLC_sum(Nmax=6)
end

#run_NLCE(ARGS, Nmax=9, Max_num_clusters=1500)
