"""
    Does the thermal average for all clusters at all the required temperatures (Temps) and fields (hs). Out put the the thermal average values of physical quantities of each order to a different cluster_info_file. Each cluster_info_file contains the thermal averages of all topologically distincit clusters in that order.

    return `multi, thermal_avg_folder_full_path`, the first entry is the multiplicities of each clusters, to use in NLC_sum function. The second entry is path to folder with thermal average data.
"""
function thermal_avg_all_clusters_xxz(; Nmax, J_xy, J_z, g, Temps, hs, simulation_folder_full_path::String, clusters_info_path::String, diag_folder_path::String, skip_exit_files=true)

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

    # make folder if doesn't exist. Don't index the folder.
    thermal_avg_folder_full_path = joinpath(simulation_folder_full_path, thermal_avg_folder_prefix)
    !ispath(thermal_avg_folder_full_path) && mkdir(thermal_avg_folder_full_path)
    #make_indexed_folder(folder_prefix=thermal_avg_folder_prefix, folder_path=simulation_folder_full_path)

    # =====================================================
    # ======= initialize holders of termal averages =======
    # =====================================================

    # Loop over the NLCE order
    for N = 2:Nmax

        # the file to hold the results (T average of clusters)
        thermal_avg_fname = "thermal_avg_order$(N)"
        thermal_avg_file = joinpath(thermal_avg_folder_full_path, thermal_avg_fname * ".jld2")

        # collection of all the clusters' hash tags of the Nth order
        cluster_hash_tags_N = cluster_hash_tags[N]

        # total number of clusters at this order
        tot_num_clusters = length(cluster_hash_tags_N)

        # sleep for 0-2 seconds to avoid parallel programs get synchronized 
        # and crush when making folder in the next line
        sleep(rand() * 2)

        # create folder to hold data of order n, if folder doesn't exist
        thermal_avg_folder_orderN = "order$(N)"
        !ispath(joinpath(thermal_avg_folder_full_path, "order$(N)")) && mkdir(joinpath(thermal_avg_folder_full_path, "order$(N)"))

        # loop over all clusters in the Nth order
        for (cluster_ind, cluster_hash_tags) in enumerate(cluster_hash_tags_N)

            """
            # =======================================
            # ======= diagonalize the cluster =======
            # =======================================

            # read in bond information
            # since bond info start to index at 0, we need to plus 1
            cluster_bonds = [[bond[1] + 1, bond[2] + 1] for bond in bonds[cluster_hash_tags]]

            # diagonalize, get eigen values
            quantities = diagonalize_cluster_xxz(N=N, sectors_info=sectors_info, bonds=cluster_bonds, J_xy=1.0, J_z=J_z / J_xy)
            """

            # =============================================
            # ======= thermally average the cluster =======
            # =============================================

            # the file to hold the thermal averages of the clusters
            thermal_avg_fname = "thermal_avg_id" * cluster_hash_tags
            thermal_avg_file = joinpath(thermal_avg_folder_full_path, thermal_avg_folder_orderN, thermal_avg_fname * ".jld2")

            # If enabled, skip if thermal average data exists
            if skip_exit_files
                if isfile(thermal_avg_file)
                    #=
                    try
                        # skip diagonalization if file exists and can be opened properly
                        h5open(thermal_avg_file, "r") do file
                        end
                        # progress checker
                        if mod(cluster_ind, 1000) == 0
                            println("finished checking #$(cluster_ind) out of total $(tot_num_clusters) clusters of order $(N)")
                        end
                        continue
                    catch e
                        # make file to replace the existing file if can't be opened properly
                        jldopen(thermal_avg_file, "w", compress=false) do file
                        end
                        saving_error_message = @sprintf "Something went wrong opening .jld2 file for cluster # %d, overwritten file " (cluster_ind)
                        saving_error_message = saving_error_message * thermal_avg_file
                        println(saving_error_message)
                    end
                    =#
                    continue
                else
                    try
                        # try making file before doing diagonalization,
                        # to prevent other programs
                        # from doing the same diagonalization simultaneously 
                        jldopen(thermal_avg_file, "a+", compress=false) do file
                        end
                    catch e
                        saving_error_message = @sprintf "Something went wrong creating .jld2 file for cluster # %d" (cluster_ind)
                        println(saving_error_message)
                    end
                end
            end

            # take thermal averages, using eigen values from diag_file
            diag_name = "diagonalized_order$(N)_id" * cluster_hash_tags
            diag_file_path = joinpath(diag_folder_path, diag_name * ".jld2")

            # print indicator of saving file
            progress_message = @sprintf "Thermalizing cluster # %d, %.2f" (cluster_ind) (100 * cluster_ind / tot_num_clusters)
            progress_message = progress_message * "% of order $(N)"
            print(progress_message * "\r")

            # read in eigen values for the cluster
            #eig_vals_data = h5open(diag_file_path, "r")
            #eig_vals = Dict("E" => read(eig_vals_data["E"]), "M" => read(eig_vals_data["M"]))
            # read in eigen values for the cluster (using HDF5)
            eig_vals = h5open(diag_file_path, "r")

            # take thermal averages of all T and h
            quantities_avg = thermal_avg_hT_loop(;
                Temps=Temps,
                hs=hs,
                J=J_xy,
                g=g,
                eig_vals=eig_vals
            )

            close(eig_vals)

            # calculate & store E, M, N (N is place holder for spin models)
            Estore = quantities_avg["E"]
            Mstore = quantities_avg["M"]
            Nstore = quantities_avg["N"]
            lnZstore = log.(quantities_avg["Z"])
            # calculate C and Chi
            # subtract <E>^2 and <M>^2 to make sure it's extensive
            Esqstore = quantities_avg["Esq"] .- quantities_avg["E"] .^ 2
            Msqstore = quantities_avg["Msq"] .- quantities_avg["M"] .^ 2

            # print indicator of saving file
            print("saving data" * "\r")
            # use the try block to prevent crush when 2 machines try to save 
            # diagonalization info of the same cluster
            try
                # save eigen values to file
                # only save E and M for each state to save space
                jldopen(thermal_avg_file, "a+", compress=true) do file
                    file["E"] = Estore
                    file["M"] = Mstore
                    file["Esq"] = Esqstore
                    file["Msq"] = Msqstore
                    file["N"] = Nstore
                    file["lnZ"] = lnZstore
                end
            catch e
                saving_error_message = @sprintf "Something wrong with saving cluster # %d, probably due to existing diagonalization data of the same cluster produced from other (parallely) running machines" (cluster_ind)
                println(saving_error_message)
            end

            # progress checker
            if mod(cluster_ind, 1000) == 0
                println("finished T averaging #$(cluster_ind) of order $(N)")
            end
        end


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
