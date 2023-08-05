using NLC
using Printf

function run_diagonalization(; J_z, J_xy=1.0)

    #================================#
    #= DEFINE SIMULATION PARAMETERS =#
    #================================#

    # g factor
    g = 1.0

    # ===================================================
    # ==== make folder to hold diagonalization data  ====
    # ===================================================

    diag_folder_prefix = @sprintf "hashtag_diagonization_data_J_z[J_xy%.4f" (J_z / J_xy)

    # diag_folder_full_path = make_indexed_folder(folder_prefix=diag_folder_prefix)

    # put data to existing folder
    diag_folder_full_path = "/nfs/home/zyt329/Research/xxz/runs/hashtag_diagonization_data_J_z[J_xy16.0000-1/"

    # ==============================
    # ==== start diagonalizing  ====
    # ==============================

    diagonalize_all_clusters_xxz(; J_xy=J_xy, J_z=J_z, Nmax=14, clusters_info_path="/nfs/home/zyt329/Research/xxz/NLC_clusters_info_JSON/triangle/", diag_folder_path=diag_folder_full_path)

end

@time run_diagonalization(J_z=16.0)
