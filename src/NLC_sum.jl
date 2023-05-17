"""
Does NLCE summation, return the raw summations.
"""
function NLC_sum(; Nmax, Max_num_clusters, J_xy, J_z, g, Temps, hs, multi, clusters_info_path::String, simulation_folder_full_path::String, thermal_avg_folder_full_path::String)

    NT = length(Temps)
    Nh = length(hs)

    # O is going to hold NLCE's partial sums (sums involving clusters
    # in a particular order order) for our four properties. These are
    # Sn in Eq. (3) of https://arxiv.org/pdf/0706.3254.pdf
    # I think it needs to hold quantities for single sites at O[:,1,:]
    O = zeros(Float64, 6, Nmax + 1, NT, Nh)

    # This array is going to hold the contribution of clusters
    # from all orders
    weights = zeros(Float64, 6, Nmax + 1, Max_num_clusters, NT, Nh)

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
        # Opening the property binary files and loading the information
        # on to arrays
        cluster_info_fname = joinpath(clusters_info_path, "NLCE_1_" * string(N) * ".txt")
        cluster_info_file = open(cluster_info_fname, "r")

        # read thermal avg data to quants_store
        thermal_avg_fname = "thermal_avg_order$(N)"
        thermal_avg_file = joinpath(thermal_avg_folder_full_path, thermal_avg_fname * ".jld")
        quants_store = load(thermal_avg_file)

        Estore = quants_store["Estore"]
        Mstore = quants_store["Mstore"]
        Esqstore = quants_store["Esqstore"]
        Msqstore = quants_store["Msqstore"]
        Nstore = quants_store["Nstore"]
        lnZstore = quants_store["lnZstore"]

        # Skiping line with order number & toplogical number
        readline(cluster_info_file)

        # Tolopogy number
        c = parse(Int64, readline(cluster_info_file)) + 1 # +1 for julia indexing

        # Going through the cluster files again to read in the
        # subcluster information
        EOF = false
        while EOF == false

            line = split(readline(cluster_info_file))

            # Skiping the bond info
            for b = 1:parse(Int64, line[1])
                readline(cluster_info_file)
            end

            # Subgraph info
            sc_Num = parse(Int64, readline(cluster_info_file))
            for sg = 1:sc_Num
                line = split(readline(cluster_info_file))
                subclusterSize = parse(Int64, line[1])
                scMult = parse(Int64, line[2])
                sb_topN = parse(Int64, line[3]) + 1 # +1 for julia indexing

                # In computing contributions from clusters, we first
                # subtract the subcluster weights, except for the single
                # site subcluster
                for i = 1:6
                    weights[i, N, c, :, :] -=
                        (weights[i, subclusterSize, sb_topN, :, :] * scMult)
                end
            end

            # We then add the properties computed in ED and subtract the
            # single site contributions. See https://arxiv.org/pdf/1207.3366.pdf
            # for more details.
            weights[1, N, c, :, :] += Estore[c, :, :] - N * singleE[:, :]
            weights[2, N, c, :, :] += Mstore[c, :, :] - N * singleM[:, :]
            weights[3, N, c, :, :] += Esqstore[c, :, :] - N * singleEsq[:, :]
            weights[4, N, c, :, :] += Msqstore[c, :, :] - N * singleMsq[:, :] 
            weights[5, N, c, :, :] += Nstore[c, :, :] - N * singleN[:, :]
            weights[6, N, c, :, :] += lnZstore[c, :, :] - N * singlelnZ[:, :]

            # We are now ready to put together the partial sums, using
            # the cluster contributions and corresponding lattice constants
            O[:, N, :, :] += multi[N, c] * weights[:, N, c, :, :]

            readline(cluster_info_file)
            if occursin("Topo", readline(cluster_info_file))
                # It's the end of the cluster_info_file
                EOF = true
                close(cluster_info_file)
                break
            else
                c += 1
            end
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

    # save O to O.jld
    save(joinpath(simulation_folder_full_path, "O" * ".jld"), "O", O)

    return O

end
