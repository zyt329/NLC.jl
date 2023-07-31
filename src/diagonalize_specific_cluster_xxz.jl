"""
    We create the Hamiltonian for all clusters and diagonalize them, get the corresponding En and measurement of other quantities for the eigen states.

    By default we set J_xy=1.0 and J_z=J_z/J_xy. Results for Different J_z could be extracted by multiplying E by J_xy.
"""
function diagonalize_specific_cluster_xxz(; J_xy::Float64, J_z::Float64, N, cluster_hash_tag, diag_folder_path::String, bonds::Dict{String,Vector{Vector{Int64}}})

    println("starting diagonalizing order " * string(N) * ", cluster id" * cluster_hash_tag)

    # Generate sector info for julia script
    sectors_info = sectors_info_gen(N=N)
    # Generate Sᶻ=0 sector info when N is even
    if iseven(N)
        m0_sectors_info = m0_sectors_info_gen(sectors_info=sectors_info, N=N)
    else
        # need to pass this variable to diagonalization function
        m0_sectors_info = nothing
    end

    # the file to hold the eigen values of the cluster
    diag_name = "diagonalized_order$(N)_id" * cluster_hash_tag
    diag_file = joinpath(diag_folder_path, diag_name * ".jld2")

    # =======================================
    # ======= diagonalize the cluster =======
    # =======================================

    # read in bond information
    # since bond info start to index at 0, we need to plus 1
    cluster_bonds = [[bond[1] + 1, bond[2] + 1] for bond in bonds[cluster_hash_tag]]

    # diagonalize, get eigen values
    quantities = diagonalize_cluster_xxz(N=N, sectors_info=sectors_info, m0_sectors_info=m0_sectors_info, bonds=cluster_bonds, J_xy=1.0, J_z=J_z / J_xy)

    # use the try block to prevent crush when 2 machines try to save 
    # diagonalization info of the same cluster
    try
        # save eigen values to file
        # only save E and M for each state to save space
        jldopen(diag_file, "w", compress=true) do file
            file["E"] = quantities[1]
            file["M"] = quantities[3]
        end
    catch e
        saving_error_message = "error overwriting " * diag_file
        println(saving_error_message)
    end

    println("finished diagonalizing cluster id" * cluster_hash_tag * " of order $(N)")

end


#diagonalize_all_clusters_xxz(5)
