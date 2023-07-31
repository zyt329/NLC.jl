"""
    Does the thermal average for all clusters at all the required temperatures (Temps) and fields (hs). Out put the the thermal average values of physical quantities of each order to a different cluster_info_file. Each cluster_info_file contains the thermal averages of all topologically distincit clusters in that order.

    return `multi, thermal_avg_folder_full_path`, the first entry is the multiplicities of each clusters, to use in NLC_sum function. The second entry is path to folder with thermal average data.
"""
function thermal_avg_cluster(; N, cluster_hash_tag, J_xy, J_z, g, Temps, hs, simulation_folder_full_path::String, diag_folder_path::String, bonds::Dict{String,Vector{Vector{Int64}}})

    NT = length(Temps)
    Nh = length(hs)

    # ===================================================#
    # ======      make folder to store data        ======#
    # ===================================================#

    # make folder to hold thermal avg data
    thermal_avg_folder_prefix = "thermal_avg_data"

    # make folder if doesn't exist. Don't index the folder.
    thermal_avg_folder_full_path = joinpath(simulation_folder_full_path, thermal_avg_folder_prefix)
    !ispath(thermal_avg_folder_full_path) && mkdir(thermal_avg_folder_full_path)

    # =====================================================
    # ======= initialize holders of termal averages =======
    # =====================================================

    # create folder to hold data of order n, if folder doesn't exist
    thermal_avg_folder_orderN = "order$(N)"
    !ispath(joinpath(thermal_avg_folder_full_path, "order$(N)")) && mkdir(joinpath(thermal_avg_folder_full_path, "order$(N)"))

    # =============================================
    # ======= thermally average the cluster =======
    # =============================================

    # the file to hold the thermal averages of the clusters
    thermal_avg_fname = "thermal_avg_id" * cluster_hash_tag
    thermal_avg_file = joinpath(thermal_avg_folder_full_path, thermal_avg_folder_orderN, thermal_avg_fname * ".jld2")

    # take thermal averages, using eigen values from diag_file
    diag_name = "diagonalized_order$(N)_id" * cluster_hash_tag
    diag_file_path = joinpath(diag_folder_path, diag_name * ".jld2")

    # print indicator of saving file
    progress_message = "Thermalizing cluster " * cluster_hash_tag * " of order $(N)"
    print(progress_message * "\r")

    # initialize variables to access data inside try block
    quantities_avg = nothing
    eig_vals = nothing
    try
        eig_vals = h5open(diag_file_path, "r")

        # take thermal averages of all T and h
        quantities_avg = thermal_avg_hT_loop(;
            Temps=Temps,
            hs=hs,
            J=J_xy,
            g=g,
            eig_vals=eig_vals
        )
    catch e
        if eig_vals != nothing # if file was opened but read with error
            # close damaged file
            close(eig_vals)
            println("closed damaged diag_file of cluster id " * cluster_hash_tag)
        end

        println("re-diagonalizing cluster id" * cluster_hash_tag)
        # diagonalize the cluster again if diag_file is damaged
        diagonalize_specific_cluster_xxz(; J_xy=J_xy, J_z=J_z, N=N, cluster_hash_tag=cluster_hash_tag, diag_folder_path=diag_folder_path, bonds=bonds)

        println("re-thermalizing cluster id" * cluster_hash_tag)

        eig_vals = h5open(diag_file_path, "r")
        # take thermal averages of all T and h
        quantities_avg = thermal_avg_hT_loop(;
            Temps=Temps,
            hs=hs,
            J=J_xy,
            g=g,
            eig_vals=eig_vals
        )
    end
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
        jldopen(thermal_avg_file, "w", compress=true) do file
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

    println("finished T averaging cluster" * cluster_hash_tag * " of order $(N)")

    return nothing

end
