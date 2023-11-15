using JLD2
using NLC
using Printf


"""
    Max_num_clusters defines maximum number of clusters in the highest order. Setting it high may result in insufficient memory.
"""
function run_NLCE(; J_z_by_J_xy, Nmax)

    #================================#
    #= DEFINE SIMULATION PARAMETERS =#
    #================================#

    # temperature values
    Temps = vcat(exp10.(range(log10(0.05), log10(25), length=41)))
    NT = length(Temps)

    # field values
    hs = vcat(range(0.0, 6.0, length=41))
    Nh = length(hs)

    # z interaction
    J_z = 1.0

    # xy interactoin
    J_xy = 1.0

    # g factor
    g = 1.0

    # ==============================
    # ====  make/find folders   ====
    # ==============================

    # path to cluster info
    clusters_info_path = "/gpfs/scratch/yzhang/Research/xxz/NLC_clusters_info_JSON/triangle_clusters_no_2_order/"

    # make simulation folders
    # simulation_folder_prefix = @sprintf "triangle_based_simulation_J_z[J_xy%.4f" (J_z / J_xy)

    simulation_folder_full_path = "/gpfs/scratch/yzhang/Research/xxz/runs/triangle_based_no2_Ising_simulation-1"

    # simulation_folder_full_path = make_indexed_folder(folder_prefix=simulation_folder_prefix, folder_path="/gpfs/scratch/yzhang/Research/xxz/runs/")

    # ==============================
    # ==== thermally average  ======
    # ==============================

    # do the thermal averages and find out multiplicities of clusters
    multi, thermal_avg_folder_full_path = thermal_avg_all_clusters_Ising(Nmax=Nmax, J_z=J_z, g=1, Temps=Temps, hs=hs, simulation_folder_full_path=simulation_folder_full_path, clusters_info_path=clusters_info_path)

    # ===========================
    # ==== do the NLC sum  ======
    # ===========================

    """
    O = NLC_sum(Nmax=Nmax, J_xy=J_xy, J_z=J_z, g=1, Temps=Temps, hs=hs, multi=multi, clusters_info_path=clusters_info_path, simulation_folder_full_path=simulation_folder_full_path, thermal_avg_folder_full_path=thermal_avg_folder_full_path)
    """

end



for J_z_by_J_xy in [1]
    # run_diagonalization(J_z=Float64(J_z))
    @time run_NLCE(J_z_by_J_xy=Float64(J_z_by_J_xy), Nmax=6)
end

#run_NLCE(ARGS, Nmax=9, Max_num_clusters=1500)
