"""
    Does the thermal average for all clusters at all the required temperatures (Temps) and fields (hs). Out put the the thermal average values of physical quantities of each order to a different cluster_info_file. Each cluster_info_file contains the thermal averages of all topologically distincit clusters in that order.

    return `multi, thermal_avg_folder_full_path`, the first entry is the multiplicities of each clusters, to use in NLC_sum function. The second entry is path to folder with thermal average data.
"""
function thermal_avg_cluster_Ising(; N, cluster_hash_tag, J_z=1.0, g, Temps, hs, simulation_folder_full_path::String, bonds::Dict{String,Vector{Vector{Int64}}})

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

    # =============================================
    # =======  Get E and M of the cluster   =======
    # =============================================

    cluster_bonds = [[bond[1] + 1, bond[2] + 1] for bond in bonds[cluster_hash_tag]]

    eig_vals = Dict("M" => Float64[], "E" => Float64[])
    for state in 0:(2^N-1)
        m = chk_m(state, N)
        push!(eig_vals["M"], (N - 2m) * 1 / 2)

        E = 0
        state_binary = digits!(zeros(Int64, 64), state, base=2)
        for bond in cluster_bonds #bond=[s1, s2], where s1,s2 are the two sites of the bond
            s1 = bond[1]
            s2 = bond[2]
            if state_binary[s1] == state_binary[s2]
                E += (1 / 4) * J_z
            else
                E += -(1 / 4) * J_z
            end
        end
        push!(eig_vals["E"], E)
    end

    # print indicator of saving file
    progress_message = "Thermalizing cluster " * cluster_hash_tag * " of order $(N)"
    print(progress_message * "\r")

    # take thermal averages of all T and h
    quantities_avg = thermal_avg_hT_loop(;
        Temps=Temps,
        hs=hs,
        J=J_z,
        g=g,
        eig_vals=eig_vals,
        Ising=true
    )

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
        saving_error_message = @sprintf "Something wrong with saving cluster # %d, probably due to existing thermal data of the same cluster produced from other (parallely) running machines" (cluster_ind)
        println(saving_error_message)
    end

    println("finished T averaging cluster" * cluster_hash_tag * " of order $(N)")

    return nothing

end
