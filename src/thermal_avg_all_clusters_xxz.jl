"""
    Does the thermal average for all clusters at all the required temperatures (Temps) and fields (hs). Out put the the thermal average values of physical quantities of each order to a different cluster_info_file. Each cluster_info_file contains the thermal averages of all topologically distincit clusters in that order.

    return `multi, thermal_avg_folder_full_path`, the first entry is the multiplicities of each clusters, to use in NLC_sum function. The second entry is path to folder with thermal average data.
"""
function thermal_avg_all_clusters_xxz(; Nmax, J_xy, J_z, g, Temps, hs, simulation_folder_full_path::String, clusters_info_path::String)

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


    # ===================================================#
    # ======      make folder to store data        ======#
    # ===================================================#

    # make folder to hold thermal avg data
    thermal_avg_folder_prefix = "thermal_avg_data"

    thermal_avg_folder_full_path = make_indexed_folder(folder_prefix=thermal_avg_folder_prefix, folder_path=simulation_folder_full_path)

    # =====================================================
    # ======= initialize holders of termal averages =======
    # =====================================================
    Estore = Dict{String,Array{Float64,2}}()
    Mstore = Dict{String,Array{Float64,2}}()
    Esqstore = Dict{String,Array{Float64,2}}()
    Msqstore = Dict{String,Array{Float64,2}}()
    Nstore = Dict{String,Array{Float64,2}}()
    lnZstore = Dict{String,Array{Float64,2}}()
    # Loop over the NLCE order
    for N = 2:Nmax

        # Generate sector info for use in diagonalization
        sectors_info = sectors_info_gen(N=N)

        # the file to hold the results (T average of clusters)
        thermal_avg_fname = "thermal_avg_order$(N)"
        thermal_avg_file = joinpath(thermal_avg_folder_full_path, thermal_avg_fname * ".jld2")

        # collection of all the clusters' hash tags of the Nth order
        cluster_hash_tags_N = cluster_hash_tags[N]

        # total number of clusters at this order
        tot_num_clusters = length(cluster_hash_tags_N)

        # loop over all clusters in the Nth order
        for (cluster_ind, cluster_hash_tags) in enumerate(cluster_hash_tags_N)

            # =======================================
            # ======= diagonalize the cluster =======
            # =======================================

            # read in bond information
            # since bond info start to index at 0, we need to plus 1
            cluster_bonds = [[bond[1] + 1, bond[2] + 1] for bond in bonds[cluster_hash_tags]]

            # diagonalize, get eigen values
            quantities = diagonalize_cluster_xxz(N=N, sectors_info=sectors_info, bonds=cluster_bonds, J_xy=1.0, J_z=J_z / J_xy)

            # =============================================
            # ======= thermally average the cluster =======
            # =============================================

            # take thermal averages, using eigen values obtained above
            quantities_avg = thermal_avg_hT_loop(;
                Temps=Temps,
                hs=hs,
                J=J_xy,
                g=g,
                quantities=quantities
            )

            # normalize by Z
            numerator_M = quantities_avg[4]
            numerator_E = quantities_avg[2]
            numerator_Esq = quantities_avg[3]
            numerator_Msq = quantities_avg[5]
            numerator_N = quantities_avg[6]
            denominator = quantities_avg[1]

            # calculate & store E, M, N (N is place holder for spin models)
            Estore[cluster_hash_tags] = numerator_E ./ denominator
            Mstore[cluster_hash_tags] = numerator_M ./ denominator
            Nstore[cluster_hash_tags] = numerator_N ./ denominator
            lnZstore[cluster_hash_tags] = log.(denominator)

            # calculate C and Chi
            # subtract <E>^2 and <M>^2 to make sure it's extensive
            Esqstore[cluster_hash_tags] = numerator_Esq ./ denominator - (numerator_E ./ denominator) .^ 2
            Msqstore[cluster_hash_tags] = numerator_Msq ./ denominator - (numerator_M ./ denominator) .^ 2

            # progress checker
            if mod(cluster_ind, 1000) == 0
                println("finished T averaging #$(cluster_ind) of order $(N)")
            end
        end

        # save T averages of clusters
        save(
            thermal_avg_file,
            "Estore",
            Estore,
            "Mstore",
            Mstore,
            "Esqstore",
            Esqstore,
            "Msqstore",
            Msqstore,
            "Nstore",
            Nstore,
            "lnZstore",
            lnZstore,
        )

        # progress checker for each order
        println("finisehd T average of ORDER ", N)
    end

    # =====================================
    # ==== save simulation parameters  ====
    # =====================================

    simulation_parameters_fname = "simulation_parameters"

    simulation_parameters_file = joinpath(simulation_folder_full_path, simulation_parameters_fname * ".jld2")

    save(simulation_parameters_file, "Nmax", Nmax,
        "Temps", Temps, "hs", hs, "J_xy", J_xy, "J_z", J_z, "k_B", k_B, "N_A", N_A, "mu_B", mu_B, "g", g)

    # create null file whose name gives parameters
    null_file_name = @sprintf "T_[%.2f-%.2f_%i]_h_[%.2f-%.2f_%i]_Kb_%.2f_g_%.4f" Temps[1] Temps[end] length(Temps) hs[1] hs[end] length(hs) k_B g

    null_file = joinpath(simulation_folder_full_path, null_file_name * ".jld2")

    save(null_file, "null_file", "null file whose name reveals parameters")

    return multi, thermal_avg_folder_full_path

end
