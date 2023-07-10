using NLC
using JLD2
using JSON

"""
    Does the thermal average for all clusters at all the required temperatures (Temps) and fields (hs). Out put the the thermal average values of physical quantities of each order to a different cluster_info_file. Each cluster_info_file contains the thermal averages of all topologically distincit clusters in that order.

    return `multi, thermal_avg_folder_full_path`, the first entry is the multiplicities of each clusters, to use in NLC_sum function. The second entry is path to folder with thermal average data.
"""
function thermal_avg_all_clusters_xxz(; Nmax, J_xy, J_z, g, Temps, hs, simulation_folder_full_path::String, clusters_info_path::String, diag_folder_path::String)

    # ===================================================#
    # ======      load up cluster information      ======#
    # ===================================================#

    # load up "lattice constants" (multiplicities) in `multi`
    # load up all the clusters' hash tags to `cluster_hash_tags`
    # load up bond information to 'bonds'
    cluster_hash_tags = []
    multi = Dict{String,Int}()
    bonds = Dict{String,Vector{Vector{Int64}}}()
    for N in 1:Nmax
        cluster_multi_info_fname = joinpath(clusters_info_path, "graph_mult_triangle_" * string(N) * ".json")

        cluster_bond_info_fname = joinpath(clusters_info_path, "graph_bond_triangle_" * string(N) * ".json")

        multi_N = JSON.parsefile(cluster_multi_info_fname)

        bonds_N = JSON.parsefile(cluster_bond_info_fname)

        merge!(multi, multi_N)
        println("Finished read in multiplicity info of order $(N)")

        merge!(bonds, bonds_N)
        println("Finished read in bond info of order $(N)")

        push!(cluster_hash_tags, keys(multi_N))
    end




    # =====================================================
    # ======= initialize holders of termal averages =======
    # =====================================================
    t = 0
    # Loop over the NLCE order
    for N = 2:Nmax

        # collection of all the clusters' hash tags of the Nth order
        cluster_hash_tags_N = cluster_hash_tags[N]

        # total number of clusters at this order
        tot_num_clusters = length(cluster_hash_tags_N)


        # loop over all clusters in the Nth order
        for (cluster_ind, cluster_hash_tags) in enumerate(cluster_hash_tags_N)
            # println("diagonalizing cluster # $cluster_ind of order $N")

            # =============================================
            # ======= thermally average the cluster =======
            # =============================================


            # take thermal averages, using eigen values from diag_file
            diag_name = "diagonalized_order$(N)_id" * cluster_hash_tags
            diag_file_path = joinpath(diag_folder_path, diag_name * ".jld2")

            dt = @elapsed begin
                eigen_vals = jldopen(diag_file_path, "r")
            end
            t += dt

            # passing in quantities
            E::Array{Float64,1} = eigen_vals["E"]
            M::Array{Float64,1} = eigen_vals["M"]
            N_tot::Array{Float64,1} = ones(size(E))

            Esq = E .^ 2
            Msq = M .^ 2


            # progress checker
            if mod(cluster_ind, 1000) == 0
                println("finished T averaging #$(cluster_ind) of order $(N)")
            end
        end



        # progress checker for each order
        println("finisehd T average of ORDER ", N)
    end


    println("Total loading time is $t")
end

clusters_info_path = "/nfs/home/zyt329/Research/xxz/NLC_clusters_info_JSON/triangle/"


@time thermal_avg_all_clusters_xxz(; Nmax=9, J_xy=1, J_z=1, g=1, Temps=nothing, hs=nothing, simulation_folder_full_path="nothing", clusters_info_path=clusters_info_path, diag_folder_path="/nfs/home/zyt329/Research/xxz/runs/hashtag_diagonization_data_J_z[J_xy1.0000-1/")