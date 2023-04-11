module NLC

# import external dependencies
using LinearAlgebra
using SparseArrays
using Arpack
using JLD
using DelimitedFiles

#=
Name changes in functions:
in xxz.jl:
    H_sector => XXZ_H_sector

in measurements_Heisbg_Hub_compare.jl:
    file name measurements_Heisbg_Hub_compare => measurements_xxz.jl
    E_Quants => diagonalize_cluster
    thermal_avg_loop => thermal_avg_hT_loop

in Diagonalization_clusters_Heisbg.jl:
    file name Diagonalization_clusters_Heisbg.jl => diagonalize_all_clusters_xxz.jl
    Diag_Heisbg => diagonalize_all_clusters_xxz

in Heisbg_Hub_compare_single_site.jl:
    file name Heisbg_Hub_compare_single_site.jl => single_site.jl
    single_site_quantities => single_site_quantities_xxz
=#


# =================================================
# ===========     Define constants    =============
# =================================================

include("consts.jl")

# =================================================
# ===========     load utilities      =============
# =================================================

include("utilities.jl")
export printing, reading_quantities

# =================================================
# ==== load functions to construct Hamiltonians ===
# =================================================

include("single_site.jl")
export single_site_quantities_Heisbg, single_site_quantities_xxz

include("xxz.jl")
export chk_m, sectors_info_gen, update_val, H_sector

include("measurements_xxz.jl")
export diagonalize_cluster, thermal_avg, thermal_avg_hT_loop

include("diagonalize_all_clusters_xxz")
export diagonalize_all_clusters_xxz













end