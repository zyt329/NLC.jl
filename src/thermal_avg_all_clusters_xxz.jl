"""
    Does the thermal average for all clusters at all the required temperatures (Temps) and fields (hs). Out put the the thermal average values of physical quantities of each order to a different cluster_info_file. Each cluster_info_file contains the thermal averages of all topologically distincit clusters in that order.

    return `multi, thermal_avg_folder_full_path`, the first entry is the multiplicities of each clusters, to use in NLC_sum function. The second entry is path to folder with thermal average data.
"""
function thermal_avg_all_clusters_xxz(; Nmax, Max_num_clusters, J_xy, J_z, g, Temps, hs, diag_folder_full_path::String, simulation_folder_full_path::String, clusters_info_path::String)

    # Holding "lattice constants" (multiplicities) for all clusters
    # in all orders. Probably good for up to order 9. Need to check
    # what is needed exactly for the size of the second dimension.
    multi = zeros(Int64, Nmax + 1, Max_num_clusters)

    #name = "Heisbg_order=" * string(12) #string(Nmax) #always read from Nmax=12th order data

    # Number of temperatures
    NT = length(Temps)

    # Number of magnetic field h
    Nh = length(hs)

    # make folder to hold thermal avg data
    thermal_avg_folder_prefix = "thermal_avg_data"

    thermal_avg_folder_full_path = make_indexed_folder(folder_prefix=thermal_avg_folder_prefix, folder_path=simulation_folder_full_path)

    # Counter for isomorphically distinct clusters
    topN = 0

    # Loop over the NLCE order
    for N = 2:Nmax

        # Initializing arrays and openning files
        Estore = zeros(Float64, Max_num_clusters, NT, Nh)
        Mstore = zeros(Float64, Max_num_clusters, NT, Nh)
        Esqstore = zeros(Float64, Max_num_clusters, NT, Nh)
        Msqstore = zeros(Float64, Max_num_clusters, NT, Nh)
        Nstore = zeros(Float64, Max_num_clusters, NT, Nh)
        lnZstore = zeros(Float64, Max_num_clusters, NT, Nh)

        # open the cluster information cluster_info_file
        # Change 1 to 2 in the line below and in the next cell to include n.n.n. bonds
        cluster_info_fname = joinpath(clusters_info_path, "NLCE_1_" * string(N) * ".txt")
        cluster_info_file = open(cluster_info_fname, "r")

        # the file to hold thermally average quantities
        thermal_avg_fname = "thermal_avg_order$(N)"
        thermal_avg_file = joinpath(thermal_avg_folder_full_path, thermal_avg_fname * ".jld")

        # Skips line with order number
        readline(cluster_info_file)

        # Going through the cluster_info_file for each order line by line
        # and reading in the cluster information
        topN = parse(Int64, readline(cluster_info_file))
        println("ORDER", N)

        EOF = false
        while EOF == false
            # NTOP starts from 1 for convenience of indexing in julia
            # topN starts from 0 as used in the cluster data from Ehsan
            NTOP = topN + 1
            line = split(readline(cluster_info_file))

            # Get the number of bonds from this line
            nB = parse(Int64, line[1])

            # Skip the lines with bond information
            for b = 1:nB
                line = split(readline(cluster_info_file))
                # next(cluster_info_file)
            end

            # Doing and storing quantieis:
            #
            # 1- Average energy, <H>
            # 2- Average magnetization, <M>
            # 3- Average energy squared, <H^2>, and
            # 4- Average magnetization squared, <M^2>
            # 5- Average total number of particle <N>

            # =============================
            # === doing thermal sums ======
            # =============================

            # reading diagonalization data
            quantities =
                reading_quantities(name="diagonalized", NTOP=topN, N=N, folder_path=diag_folder_full_path)

            # calculate thermal averages for all h and T
            quantities_avg = thermal_avg_hT_loop(;
                Temps=Temps,
                hs=hs,
                J=J_xy,
                g=g,
                quantities=quantities
            )

            numerator_M = quantities_avg[4]
            numerator_E = quantities_avg[2]
            numerator_Esq = quantities_avg[3]
            numerator_Msq = quantities_avg[5]
            numerator_N = quantities_avg[6]
            denominator = quantities_avg[1]

            # end = time.time()
            # print(end - start)
            Estore[NTOP, :, :] = numerator_E ./ denominator
            Mstore[NTOP, :, :] = numerator_M ./ denominator
            Nstore[NTOP, :, :] = numerator_N ./ denominator
            lnZstore[NTOP, :, :] = log.(denominator)

            # It is important to do the following subtractions at the
            # individual cluster level to make the quantities extensive
            # for the NLCE.
            Esqstore[NTOP, :, :] = numerator_Esq ./ denominator - Estore[NTOP, :, :] .^ 2
            Msqstore[NTOP, :, :] = numerator_Msq ./ denominator - Mstore[NTOP, :, :] .^ 2

            # ---------------------------------------------------------

            # Here, we take the opportunity to read in the "lattice constants"
            # (multiplicities) for each topological cluster if we have reached
            # the end of the cluster_info_file.

            # skipping the subgraph information for now
            # print(line)
            sc_Num = parse(Int64, readline(cluster_info_file))
            for s = 1:sc_Num
                readline(cluster_info_file)
            end
            readline(cluster_info_file)

            # Checking if we have reached the end of the cluster_info_file
            if occursin("Topo", readline(cluster_info_file))
                # next(cluster_info_file)
                for i = 1:(topN+1)
                    multi[N, i] = parse(Int64, split(readline(cluster_info_file))[2])
                end
                EOF = true
                close(cluster_info_file)
                break
            else
                # If not yet done with clusters
                topN += 1
                if mod(topN, 1000) == 0
                    println("finished T averaging topN=$(topN)")
                end
            end
        end
        # Saving the properties to files
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
    end

    # =====================================
    # ==== save simulation parameters  ====
    # =====================================

    simulation_parameters_fname = "simulation_parameters"

    simulation_parameters_file = joinpath(simulation_folder_full_path, simulation_parameters_fname * ".jld")

    save(simulation_parameters_file, "Nmax", Nmax,
        "Temps", Temps, "hs", hs, "J_xy", J_xy, "J_z", J_z, "k_B", k_B, "N_A", N_A, "mu_B", mu_B, "g", g)

    return multi, thermal_avg_folder_full_path

end
