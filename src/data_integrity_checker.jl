"""
    We create the Hamiltonian for all clusters and diagonalize them, get the corresponding En and measurement of other quantities for the eigen states.

    By default we set J_xy=1.0 and J_z=J_z/J_xy. Results for Different J_z could be extracted by multiplying E by J_xy.
"""
function data_integrity_checker(;Nmax::Int64, clusters_info_path::String, diag_folder_path::String)

    # ===================================================#
    # ======      load up cluster information      ======#
    # ===================================================#

    # load up all the clusters' hash tags to `cluster_hash_tag`
    # load up bond information to 'bonds'
    cluster_hash_tags = []
    bonds = Dict{String,Vector{Vector{Int64}}}()
    for N in 1:Nmax
        cluster_bond_info_fname = joinpath(clusters_info_path, "graph_bond_triangle_" * string(N) * ".json")

        bonds_N = JSON.parsefile(cluster_bond_info_fname)

        merge!(bonds, bonds_N)
        println("Finished read in bond info of order $(N)")

        push!(cluster_hash_tags, keys(bonds_N))
    end

    # prefix of files to hold eigenvalues
    name = "xxz_eigvals"

    # Loop over the NLCE order
    for N in 2:Nmax

        println("starting checking order " * string(N))

        """
        # Generate sector info for julia script
        sectors_info = sectors_info_gen(N=N)
        # Generate Sᶻ=0 sector info when N is even
        if iseven(N)
            m0_sectors_info = m0_sectors_info_gen(sectors_info=sectors_info, N=N)
        else
            # need to pass this variable to diagonalization function
            m0_sectors_info = nothing
        end"""

        # collection of all the clusters' hash tags of the Nth order
        cluster_hash_tags_N = cluster_hash_tags[N]

        # total number of clusters at this order
        tot_num_clusters = length(cluster_hash_tags_N)

        # loop over all clusters in the Nth order
        for (cluster_ind, cluster_hash_tag) in enumerate(cluster_hash_tags_N)

            # the file to hold the eigen values of the cluster
            diag_name = "diagonalized_order$(N)_id" * cluster_hash_tag
            diag_file = joinpath(diag_folder_path, diag_name * ".jld2")
          
            # print indicator of checking file
            progress_message =  "Checking cluster ID " * cluster_hash_tag * @sprintf ", # %d, %.2f" (cluster_ind) (100 * cluster_ind / tot_num_clusters)
            progress_message = progress_message * "% of order $(N)"
            print(progress_message * "\r")

            # if file exists
            # delete the file if it's compromised
            if isfile(diag_file)
                try
                    # try making file before doing diagonalization,
                    # to prevent other programs
                    # from doing the same diagonalization simultaneously 
                    data = load(diag_file) 

                    # dimension of Hilbert space
                    H_dim = 2^N

                    # check number of eigenvalues
                    @assert(length(data["E"])==H_dim, "Not enough eigenvalues stored for E")

                    # check number of eigenvalues
                    @assert(length(data["M"])==H_dim, "Not enough eigenvalues stored for M")
                catch e
                    # if something wrong with the data
                    # try opening with other method
                    try
                        data = jldopen(diag_file, "r")

                        # dimension of Hilbert space
                        H_dim = 2^N

                        # check number of eigenvalues
                        @assert(length(data["E"])==H_dim, "Not enough eigenvalues stored for E")

                        # check number of eigenvalues
                        @assert(length(data["M"])==H_dim, "Not enough eigenvalues stored for M")
                    catch e
                        #@warn  e

                        # delete the compormised file
                        rm(diag_file)

                        error_message = "Deleted: Incomplete data for cluster # " * cluster_hash_tag
                        println(error_message)
                    end
                end
            end
            

            """
            # =======================================
            # ======= diagonalize the cluster =======
            # =======================================

            # read in bond information
            # since bond info start to index at 0, we need to plus 1
            cluster_bonds = [[bond[1] + 1, bond[2] + 1] for bond in bonds[cluster_hash_tag]]

            # print indicator of saving file
            progress_message = @sprintf "diagonalizing cluster # %d, %.2f" (cluster_ind) (100 * cluster_ind / tot_num_clusters)
            progress_message = progress_message * "% of order $(N)"
            print(progress_message * "\r")

            # diagonalize, get eigen values
            quantities = diagonalize_cluster_xxz(N=N, sectors_info=sectors_info, m0_sectors_info=m0_sectors_info, bonds=cluster_bonds, J_xy=1.0, J_z=J_z / J_xy)

            # print indicator of saving file
            print("saving data" * "\r")
            # use the try block to prevent crush when 2 machines try to save 
            # diagonalization info of the same cluster
            try
                # save eigen values to file
                # only save E and M for each state to save space
                jldopen(diag_file, "a+", compress=true) do file
                    file["E"] = quantities[1]
                    file["M"] = quantities[3]
                end
            catch e
                saving_error_message = @sprintf "Something wrong with saving cluster # %d, probably due to existing diagonalization data of the same cluster produced from other (parallely) running machines" (cluster_ind)
                println(saving_error_message)
            end
            """



            # progress checker
            if mod(cluster_ind, 1000) == 0
                println("finished checking #$(cluster_ind) of order $(N)")
            end

        end

        println("finished checking data of clusters of order " * string(N))
    end
end


#diagonalize_all_clusters_xxz(5)
