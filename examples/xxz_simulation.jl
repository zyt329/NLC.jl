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
    Temps = vcat(exp10.(range(log10(0.05), log10(5), length=41)))
    NT = length(Temps)

    # field values
    hs = vcat(range(0.0, 6.0, length=41))
    Nh = length(hs)

    # Jz/Jxy ratio
    J_z_by_J_xy = J_z_by_J_xy

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
    clusters_info_path = "/nfs/home/zyt329/Research/xxz/NLC_clusters_info_JSON/triangle/"

    # make simulation folders
    # simulation_folder_prefix = @sprintf "hashtag_simulation_J_z[J_xy%.4f" (J_z / J_xy)

    simulation_folder_full_path = "/nfs/home/zyt329/Research/xxz/runs/hashtag_simulation_J_z[J_xy32.0000-2"

    # simulation_folder_full_path = make_indexed_folder(folder_prefix=simulation_folder_prefix, folder_path="/nfs/home/zyt329/Research/xxz/runs/")

    # ==============================
    # ==== thermally average  ======
    # ==============================

    # do the thermal averages and find out multiplicities of clusters
    multi, thermal_avg_folder_full_path = thermal_avg_all_clusters_xxz(Nmax=Nmax, J_xy=J_xy, J_z=J_z, g=1, Temps=Temps, hs=hs, simulation_folder_full_path=simulation_folder_full_path, clusters_info_path=clusters_info_path, diag_folder_path="/nfs/home/zyt329/Research/xxz/runs/hashtag_diagonization_data_J_z[J_xy32.0000-1/")

    # ===========================
    # ==== do the NLC sum  ======
    # ===========================

    """
    O = NLC_sum(Nmax=Nmax, J_xy=J_xy, J_z=J_z, g=1, Temps=Temps, hs=hs, multi=multi, clusters_info_path=clusters_info_path, simulation_folder_full_path=simulation_folder_full_path, thermal_avg_folder_full_path=thermal_avg_folder_full_path)
    """

end



for J_z_by_J_xy in [32]
    # run_diagonalization(J_z=Float64(J_z))
    @time run_NLCE(J_z_by_J_xy=Float64(J_z_by_J_xy), Nmax=14)
end

#run_NLCE(ARGS, Nmax=9, Max_num_clusters=1500)
