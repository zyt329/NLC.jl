"""
Does NLCE summation, return the raw summations.
"""
function NLC_sum_Ising(; Nmax, clusters_info_path::String, simulation_folder_full_path::String, thermal_avg_folder_full_path::String, print_progress=false)

    simulation_param_file_path = joinpath(simulation_folder_full_path, "simulation_parameters.jld2")
    simulation_param_file = h5open(simulation_param_file_path, "r")

    Temps = read(simulation_param_file["Temps"])
    hs = read(simulation_param_file["hs"])
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

    # To hold the weights of all clusters up to the order of (Nmax - 1)
    weights = Dict()
    for N in 1:(Nmax-1)
        for cluster_hash_tag in cluster_hash_tags[N]
            weights[cluster_hash_tag] = zeros(Float64, 6, NT, Nh)
        end
    end

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
            if print_progress
                print(progress_message * "\r")
            end

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

                thermal_avg_cluster_Ising(N=N, cluster_hash_tag=cluster_hash_tag, J_z=1.0, g=g, Temps=Temps, hs=hs, simulation_folder_full_path=simulation_folder_full_path, bonds=bonds)

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

            # initialize weights for the current cluster if it's the last order 
            # This saves memory by not storing all weights of the last order 
            if N == Nmax
                weights_Nmax = zeros(Float64, 6, NT, Nh)
            end

            # loop over sub clusters of the current cluster
            for sub_cluster_hash_tag in keys(subgraph_multi[cluster_hash_tag])

                scMult = subgraph_multi[cluster_hash_tag][sub_cluster_hash_tag]

                if N != Nmax # for all but the last order:
                    weights[cluster_hash_tag][:, :, :] -=
                        (weights[sub_cluster_hash_tag][:, :, :] * scMult)
                else # for the last order:
                    weights_Nmax[:, :, :] -=
                        (weights[sub_cluster_hash_tag][:, :, :] * scMult)
                end
            end
            
            # check the number of sites (necessary for triangle-based expansion)
            cluster_bonds = [[bond[1] + 1, bond[2] + 1] for bond in bonds[cluster_hash_tag]]
            N_sites = num_sites(cluster_bonds)

            # subtract single-site contributions
            if N != Nmax # for all but the last order:
                weights[cluster_hash_tag][1, :, :] += Estore[:, :] - N_sites * singleE[:, :]
                weights[cluster_hash_tag][2, :, :] += Mstore[:, :] - N_sites * singleM[:, :]
                weights[cluster_hash_tag][3, :, :] += Esqstore[:, :] - N_sites * singleEsq[:, :]
                weights[cluster_hash_tag][4, :, :] += Msqstore[:, :] - N_sites * singleMsq[:, :]
                weights[cluster_hash_tag][5, :, :] += Nstore[:, :] - N_sites * singleN[:, :]
                weights[cluster_hash_tag][6, :, :] += lnZstore[:, :] - N_sites * singlelnZ[:, :]
            else # for the last order:
                weights_Nmax[1, :, :] += Estore[:, :] - N_sites * singleE[:, :]
                weights_Nmax[2, :, :] += Mstore[:, :] - N_sites * singleM[:, :]
                weights_Nmax[3, :, :] += Esqstore[:, :] - N_sites * singleEsq[:, :]
                weights_Nmax[4, :, :] += Msqstore[:, :] - N_sites * singleMsq[:, :]
                weights_Nmax[5, :, :] += Nstore[:, :] - N_sites * singleN[:, :]
                weights_Nmax[6, :, :] += lnZstore[:, :] - N_sites * singlelnZ[:, :]
            end

            # We are now ready to put together the partial sums, using
            # the cluster contributions and corresponding lattice constants
            if N != Nmax # for all but the last order:
                O[:, N, :, :] += multi[cluster_hash_tag] * weights[cluster_hash_tag][:, :, :]
            else # for the last order:
                O[:, N, :, :] += multi[cluster_hash_tag] * weights_Nmax[:, :, :]
            end

            # progress checker
            if mod(cluster_ind, 1000) == 0
                println("finished NLC sum of #$(cluster_ind) of order $(N)")
            end
        end

        println("\n finished NLCE order" * string(N))
    end

    # save O to O.jld2
    save(joinpath(simulation_folder_full_path, "O" * ".jld2"), "O", O)

    return O

end
