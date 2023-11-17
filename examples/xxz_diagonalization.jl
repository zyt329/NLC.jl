using NLC
using Printf

function run_diagonalization(; J_z, J_xy=1.0, Nmax)

    #================================#
    #= DEFINE SIMULATION PARAMETERS =#
    #================================#

    # g factor
    g = 1.0

    # ===================================================
    # ==== make folder to hold diagonalization data  ====
    # ===================================================

    # set up the prefix of the folder to hold the diagonalization data
    diag_folder_prefix = @sprintf "test_diagonization_data_J_z[J_xy%.4f" (J_z / J_xy)

    # make folder for diagonalization data automatically
    diag_folder_full_path = make_indexed_folder(folder_prefix=diag_folder_prefix)

    # ==============================
    # ==== start diagonalizing  ====
    # ==============================

    # input cluster info path
    clusters_info_path = "../cluster_info/triangle/"

    diagonalize_all_clusters_xxz(; J_xy=J_xy, J_z=J_z, Nmax=Nmax, clusters_info_path=clusters_info_path, diag_folder_path=diag_folder_full_path)

end

# assuming J_xy = 1.0

# J_z = 1.0 would be Heisenberg model
@time run_diagonalization(J_z=1.0, Nmax=9)
