"""
Does NLCE summation, return the raw summations.
"""
function NLC_sum(; Nmax, clusters_info_path::String, simulation_folder_full_path::String, thermal_avg_folder_full_path::String, diag_folder_path::String)

    simulation_param_file_path = joinpath(simulation_folder_full_path, "simulation_parameters.jld2")
    simulation_param_file = h5open(simulation_param_file_path, "r")

    Temps = read(simulation_param_file["Temps"])
    hs = read(simulation_param_file["hs"])
    J_xy = read(simulation_param_file["J_xy"])
    J_z = read(simulation_param_file["J_z"])
    g = read(simulation_param_file["g"])

    NT = length(Temps)
    Nh = length(hs)

    close(simulation_param_file)

    # ===================================================#
    # ======      load up cluster information      ======#
    # ===================================================#

    # load up "lattice constants" (multiplicities) in `multi`
    # load up all the clusters' hash tags to `cluster_hash_tags`
    # load up bond information to 'bonds'
    cluster_hash_tags = []
    multi = Dict{String,Int}()
    subgraph_multi = Dict{String,Any}()
    bonds = Dict{String,Vector{Vector{Int64}}}()
    for N in 1:Nmax
        cluster_multi_info_fname = joinpath(clusters_info_path, "graph_mult_triangle_" * string(N) * ".json")

        subcluster_multi_info_fname = joinpath(clusters_info_path, "subgraph_mult_triangle_" * string(N) * ".json")

        cluster_bond_info_fname = joinpath(clusters_info_path, "graph_bond_triangle_" * string(N) * ".json")

        multi_N = JSON.parsefile(cluster_multi_info_fname)

        bonds_N = JSON.parsefile(cluster_bond_info_fname)

        merge!(bonds, bonds_N)
        println("Finished read in bond info of order $(N)")

        subgraph_multi_N = JSON.parsefile(subcluster_multi_info_fname)

        merge!(multi, multi_N)
        println("Finished read in multiplicity info of order $(N)")

        merge!(subgraph_multi, subgraph_multi_N)
        println("Finished read in subgraph info of order $(N)")

        push!(cluster_hash_tags, keys(multi_N))
    end

    #= debug
    println("Size of cluster_hash_tags is $(Base.summarysize(cluster_hash_tags)/10^6)")
    println("Size of multi is $(Base.summarysize(multi)/10^6)")
    println("Size of subgraph_multi is $(Base.summarysize(subgraph_multi)/10^6)")
    =#

    # get cluster_hash_tag of single site
    single_site_hash_tag = cluster_hash_tags[1]

    # get all cluster_hash_tags
    all_cluster_hash_tags = keys(multi)

    # =======================================================#
    # ======      Initialize holder of quantities      ======#
    # =======================================================#

    # O is going to hold NLCE's partial sums (sums involving clusters
    # in a particular order order) for our four properties. These are
    # Sn in Eq. (3) of https://arxiv.org/pdf/0706.3254.pdf
    # I think it needs to hold quantities for single sites at O[:,1,:]
    O = zeros(Float64, 6, Nmax + 1, NT, Nh)

    # To hold the weights of all clusters
    weights = Dict(cluster_hash_tag => zeros(Float64, 6, NT, Nh) for cluster_hash_tag in all_cluster_hash_tags)

    # Hard coding the contributions from the single site by hand
    # Be careful with the order of quantities in O and in julia code,
    # they might not be the same and caution is required.
    single_site_quants = single_site_quantities_xxz(Ts=Temps, hs=hs, g=g)
    singleE = single_site_quants["E"]
    O[1, 1, :, :] = singleE
    singleEsq = single_site_quants["Esq"] - single_site_quants["E"] .^ 2
    O[3, 1, :, :] = singleEsq
    singleM = single_site_quants["M"]
    O[2, 1, :, :] = singleM
    singleMsq = single_site_quants["Msq"] - single_site_quants["M"] .^ 2
    O[4, 1, :, :] = singleMsq
    singleN = single_site_quants["N"]
    O[5, 1, :, :] = singleN
    singlelnZ = log.(single_site_quants["Z"])
    O[6, 1, :, :] = singlelnZ

    # Loop over NLCE orders
    for N = 2:Nmax

        # folder name to hold order N data
        thermal_avg_folder_orderN = "order$(N)"

        # total number of clusters at this order
        tot_num_clusters_N = length(cluster_hash_tags[N])

        # loop over all clusters
        for (cluster_ind, cluster_hash_tag) in enumerate(cluster_hash_tags[N])

            # print indicator of NLC
            progress_message = @sprintf "Doing NLC sum for cluster # %d, %.2f" (cluster_ind) (100 * cluster_ind / tot_num_clusters_N)
            progress_message = progress_message * "% of order $(N)"
            print(progress_message * "\r")

            # the file holding the thermal averages of the clusters

            thermal_avg_fname = "thermal_avg_id" * cluster_hash_tag
            thermal_avg_file = joinpath(thermal_avg_folder_full_path, thermal_avg_folder_orderN, thermal_avg_fname * ".jld2")

            # initialize variables (to export it from the try block)
            Estore = nothing
            Mstore = nothing
            Esqstore = nothing
            Msqstore = nothing
            Nstore = nothing
            lnZstore = nothing
            quants_store = nothing
            try # in case thermal average file is damaged
                quants_store = h5open(thermal_avg_file, "r")
                #quants_store = jldopen(thermal_avg_file)

                Estore = read(quants_store["E"])
                Mstore = read(quants_store["M"])
                Esqstore = read(quants_store["Esq"])
                Msqstore = read(quants_store["Msq"])
                Nstore = read(quants_store["N"])
                lnZstore = read(quants_store["lnZ"])

                # close thermal_avg_file (release memory hopefully)
                close(quants_store)
            catch e # thermal average the cluster manually

                if quants_store != nothing # if file was opened but read with error
                    # close damaged file
                    close(quants_store)
                    println("closed damaged thermal_avg_file of cluster id " * cluster_hash_tag)
                end

                # print indicator of re-thermalizing cluster
                println("error reading cluster " * cluster_hash_tag * ", re-thermalizing it.")

                thermal_avg_cluster(N=N, cluster_hash_tag=cluster_hash_tag, J_xy=J_xy, J_z=J_z, g=g, Temps=Temps, hs=hs, simulation_folder_full_path=simulation_folder_full_path, diag_folder_path=diag_folder_path, bonds=bonds)

                quants_store = h5open(thermal_avg_file, "r")
                #quants_store = jldopen(thermal_avg_file)

                Estore = read(quants_store["E"])
                Mstore = read(quants_store["M"])
                Esqstore = read(quants_store["Esq"])
                Msqstore = read(quants_store["Msq"])
                Nstore = read(quants_store["N"])
                lnZstore = read(quants_store["lnZ"])

                # close thermal_avg_file (release memory hopefully)
                close(quants_store)
            end



            # initialize weights for the cluster
            # weights = zeros(Float64, 6, NT, Nh)

            # loop over sub clusters of the current cluster
            for sub_cluster_hash_tag in keys(subgraph_multi[cluster_hash_tag])

                scMult = subgraph_multi[cluster_hash_tag][sub_cluster_hash_tag]

                # In computing contributions from clusters, we first
                # subtract the subcluster weights, except for the single
                # site subcluster
                for i = 1:6 # for all 6 quantities

                    weights[cluster_hash_tag][i, :, :] -=
                        (weights[sub_cluster_hash_tag][i, :, :] * scMult)

                    #=
                    weights[i, :, :] -=
                        (weights[i, :, :] * scMult)
                        =#
                end
            end

            weights[cluster_hash_tag][1, :, :] += Estore[:, :] - N * singleE[:, :]
            weights[cluster_hash_tag][2, :, :] += Mstore[:, :] - N * singleM[:, :]
            weights[cluster_hash_tag][3, :, :] += Esqstore[:, :] - N * singleEsq[:, :]
            weights[cluster_hash_tag][4, :, :] += Msqstore[:, :] - N * singleMsq[:, :]
            weights[cluster_hash_tag][5, :, :] += Nstore[:, :] - N * singleN[:, :]
            weights[cluster_hash_tag][6, :, :] += lnZstore[:, :] - N * singlelnZ[:, :]

            #=
            weights[1, :, :] += Estore[:, :] - N * singleE[:, :]
            weights[2, :, :] += Mstore[:, :] - N * singleM[:, :]
            weights[3, :, :] += Esqstore[:, :] - N * singleEsq[:, :]
            weights[4, :, :] += Msqstore[:, :] - N * singleMsq[:, :]
            weights[5, :, :] += Nstore[:, :] - N * singleN[:, :]
            weights[6, :, :] += lnZstore[:, :] - N * singlelnZ[:, :]
            =#

            # We are now ready to put together the partial sums, using
            # the cluster contributions and corresponding lattice constants
            O[:, N, :, :] += multi[cluster_hash_tag] * weights[cluster_hash_tag][:, :, :]
            # O[:, N, :, :] += multi[cluster_hash_tag] * weights[:, :, :]

        end

        println("\n finished NLCE order" * string(N))
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
