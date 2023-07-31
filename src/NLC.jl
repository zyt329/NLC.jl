module NLC

# import external dependencies
using LinearAlgebra
using SparseArrays
using Arpack
using JLD2
using DelimitedFiles
using Dates
using Statistics
using Random
using Distributions
using OffsetArrays: Origin
using OffsetArrays
using Printf
using JSON
using HDF5


#=
Name changes in functions:
in xxz.jl:
    H_sector => XXZ_H_sector

in measurements_Heisbg_Hub_compare.jl:
    file name measurements_Heisbg_Hub_compare => measurements_xxz.jl
    E_Quants => diagonalize_cluster_xxz
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
export printing, reading_quantities, make_indexed_folder

# =================================================
# ==== load functions to construct Hamiltonians ===
# =================================================

include("xxz.jl")
export chk_m, sectors_info_gen, update_val, H_sector

# =================================================
# == load functions to diagonalize Hamiltonians  ==
# =================================================

include("measurements_xxz.jl")
export diagonalize_cluster_xxz, thermal_avg, thermal_avg_hT_loop

include("diagonalize_all_clusters_xxz.jl")
export diagonalize_all_clusters_xxz

include("diagonalize_specific_cluster_xxz.jl")
export diagonalize_specific_cluster_xxz

# =================================================
# == load functions for single-site calculation  ==
# =================================================

include("single_site.jl")
export single_site_quantities_Heisbg, single_site_quantities_xxz

# =================================================
# ====== load functions for NLC calculation  ======
# =================================================

include("thermal_avg_all_clusters_xxz.jl")
export thermal_avg_all_clusters_xxz

# thermal average single cluster
include("thermal_avg_specific_cluster_xxz.jl")
export thermal_avg_cluster

include("NLC_sum.jl")
export NLC_sum

# ===============================================
# ====== load functions to do resummation  ======
# ===============================================

include("resummation.jl")
export resummation

# =====================================================
# ====== load functions to check data integrity  ======
# =====================================================

include("data_integrity_checker.jl")
export data_integrity_checker

end
