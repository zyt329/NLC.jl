"""
Does NLCE summation, return the raw summations.
"""
function NLC_sum(; Nmax, J_xy, J_z, g, Temps, hs, multi, clusters_info_path::String, simulation_folder_full_path::String, thermal_avg_folder_full_path::String)

    NT = length(Temps)
    Nh = length(hs)

    # ===================================================#
    # ======      load up cluster information      ======#
    # ===================================================#

    # load up "lattice constants" (multiplicities) in `multi`
    # load up all the clusters' hash tags to `cluster_hash_tags`
    # load up bond information to 'bonds'
    cluster_hash_tags = []
    multi = Dict{String,Int}()
    subgraph_multi = Dict{String,Any}()
    for N in 1:Nmax
        cluster_multi_info_fname = joinpath(clusters_info_path, "graph_mult_triangle_" * string(N) * ".json")

        subcluster_multi_info_fname = joinpath(clusters_info_path, "subgraph_mult_triangle_" * string(N) * ".json")

        multi_N = JSON.parsefile(cluster_multi_info_fname)

        subgraph_multi_N = JSON.parsefile(subcluster_multi_info_fname)

        merge!(multi, multi_N)
        println("Finished read in multiplicity info of order $(N)")

        merge!(subgraph_multi, subgraph_multi_N)
        println("Finished read in bond info of order $(N)")

        push!(cluster_hash_tags, keys(multi_N))
    end

    # get cluster_hash_tag of single site
    single_site_hash_tag = cluster_hash_tags[1]

    # get all cluster_hash_tags
    all_cluster_hash_tags = keys(multi)
    # ========= construction in progress !!! ================

    # O is going to hold NLCE's partial sums (sums involving clusters
    # in a particular order order) for our four properties. These are
    # Sn in Eq. (3) of https://arxiv.org/pdf/0706.3254.pdf
    # I think it needs to hold quantities for single sites at O[:,1,:]
    O = zeros(Float64, 6, Nmax + 1, NT, Nh)

    # This array is going to hold the contribution of clusters
    # from all orders
    weights = Dict(cluster_hash_tag => zeros(Float64, 6, NT, Nh) for cluster_hash_tag in all_cluster_hash_tags)

    # Hard coding the contributions from the single site by hand
    # Be careful with the order of quantities in O and in julia code,
    # they might not be the same and caution is required.
    single_site_quants = single_site_quantities_xxz(Ts=Temps, hs=hs, g=g)
    singleE = single_site_quants[2]
    O[1, 1, :, :] = singleE
    singleEsq = single_site_quants[3] - single_site_quants[2] .^ 2
    O[3, 1, :, :] = singleEsq
    singleM = single_site_quants[4]
    O[2, 1, :, :] = singleM
    singleMsq = single_site_quants[5] - single_site_quants[4] .^ 2
    O[4, 1, :, :] = singleMsq
    singleN = single_site_quants[6]
    O[5, 1, :, :] = singleN
    singlelnZ = log.(single_site_quants[1])
    O[6, 1, :, :] = singlelnZ

    # Loop over NLCE orders
    for N = 2:Nmax
        """
        # read thermal avg data to quants_store
        thermal_avg_fname = "thermal_avg_order$(N)"
        thermal_avg_file = joinpath(thermal_avg_folder_full_path, thermal_avg_fname * ".jld2")
        quants_store = load(thermal_avg_file)

        Estore = quants_store["Estore"]
        Mstore = quants_store["Mstore"]
        Esqstore = quants_store["Esqstore"]
        Msqstore = quants_store["Msqstore"]
        Nstore = quants_store["Nstore"]
        lnZstore = quants_store["lnZstore"]
        """

        thermal_avg_folder_orderN = "order$(N)"

        # loop over all clusters
        for cluster_hash_tag in cluster_hash_tags[N]
            # Subtract sub cluster contribution
            #sc_Num = parse(Int64, readline(cluster_info_file))

            # the file holding the thermal averages of the clusters

            thermal_avg_fname = "thermal_avg_id" * cluster_hash_tag
            thermal_avg_file = joinpath(thermal_avg_folder_full_path, thermal_avg_folder_orderN, thermal_avg_fname * ".jld2")
            quants_store = jldopen(thermal_avg_file)

            Estore = quants_store["E"]
            Mstore = quants_store["M"]
            Esqstore = quants_store["Esq"]
            Msqstore = quants_store["Msq"]
            Nstore = quants_store["N"]
            lnZstore = quants_store["lnZ"]

            # loop over sub clusters of the current cluster
            for sub_cluster_hash_tag in keys(subgraph_multi[cluster_hash_tag])

                scMult = subgraph_multi[cluster_hash_tag][sub_cluster_hash_tag]

                # In computing contributions from clusters, we first
                # subtract the subcluster weights, except for the single
                # site subcluster
                for i = 1:6 # for all 6 quantities
                    weights[cluster_hash_tag][i, :, :] -=
                        (weights[sub_cluster_hash_tag][i, :, :] * scMult)
                end
            end

            weights[cluster_hash_tag][1, :, :] += Estore[:, :] - N * singleE[:, :]
            weights[cluster_hash_tag][2, :, :] += Mstore[:, :] - N * singleM[:, :]
            weights[cluster_hash_tag][3, :, :] += Esqstore[:, :] - N * singleEsq[:, :]
            weights[cluster_hash_tag][4, :, :] += Msqstore[:, :] - N * singleMsq[:, :]
            weights[cluster_hash_tag][5, :, :] += Nstore[:, :] - N * singleN[:, :]
            weights[cluster_hash_tag][6, :, :] += lnZstore[:, :] - N * singlelnZ[:, :]

            # We are now ready to put together the partial sums, using
            # the cluster contributions and corresponding lattice constants
            O[:, N, :, :] += multi[cluster_hash_tag] * weights[cluster_hash_tag][:, :, :]


        end

        println("finished NLCE order" * string(N))
    end
    # printing out partial sums for order 1 and 2, at NT=Ti
    #=
    Ti = NT
    hi = Nh
    for N in 1:Nmax
        print("E of order " * string(N) * " is " * string(O[1, N, Ti, hi]))
        print("M of order " * string(N) * " is " * string(O[2, N, Ti, hi]))
        print("Esq of order " * string(N) * " is " * string(O[3, N, Ti, hi]))
        print("Msq of order " * string(N) * " is " * string(O[4, N, Ti, hi]))
        print("N of order " * string(N) * " is " * string(O[5, N, Ti, hi]))
        print("lnZ of order " * string(N) * " is " * string(O[6, N, Ti, hi]))
        print("\n")
    end=#

    # save O to O.jld2
    save(joinpath(simulation_folder_full_path, "O" * ".jld2"), "O", O)

    return O

end
