"""
    We create the Hamiltonian for all clusters and diagonalize them, get the corresponding En and measurement of other quantities for the eigen states.

    By default we set J_xy=1.0 and J_z=J_z/J_xy. Results for Different J_z could be extracted by multiplying E by J_xy.
"""
function diagonalize_all_clusters_xxz(; J_xy::Float64, J_z::Float64, Nmax::Int64, clusters_info_path::String, diag_folder_path::String)
    # Maximum NLCE order and
    # J1 and J2 values are input arguments

    # Holding "lattice constants" (multiplicities) for all clusters
    # in all orders. Probably good for up to order 9. Need to check
    # what is needed exactly for the size of the second dimension.
    multi = zeros(Int64, Nmax + 1, 16000)

    # Holding the two site numbers and optionally the bond type
    # for every bond.
    site1 = zeros(Int64, 40)
    site2 = zeros(Int64, 40)
    typ = zeros(Int64, 40)

    # Counter for isomorphically distinct clusters
    topN = 0

    # prefix of files to hold eigenvalues
    name = "xxz_eigvals"

    # Loop over the NLCE order
    for N in 2:Nmax

        println("starting order " * string(N))

        # Generate sector info for julia script
        sectors_info = sectors_info_gen(N=N)

        # Change 1 to 2 in the line below and in the next cell to include n.n.n. bonds
        fname = joinpath(clusters_info_path, "NLCE_1_" * string(N) * ".txt")
        file = open(fname, "r")

        # name of the file storing diagonalization result
        diag_name = "diagonalized"

        # Skips line with order number
        readline(file)

        # Going through the file for each order line by line and reading in the cluster information
        topN = parse(Int64, readline(file))

        EOF = false
        while EOF == false

            # passing topN to julia script
            NTOP = topN + 1
            if mod(topN, 1000) == 0
                println("diagonalizing order" * string(N) * "topo#" * string(topN))
            end

            line = split(readline(file))

            # Get the number of bonds from this line, and pass it to julia script
            nB = parse(Int64, line[1])

            # Read in the bond information, and pass it to julia script
            bonds = Array{Int64,1}[]
            for b in 1:nB
                line = split(readline(file))
                # typ[b] = parse(Int64, line[3])
                site1[b] = parse(Int64, line[1])
                site2[b] = parse(Int64, line[2])
                push!(bonds, [site1[b] + 1, site2[b] + 1])
            end

            # calculate eigenvalues of matrix and expectations of quantities in eigenstates. By default we set J_xy=1.0 and J_z=J_z/J_xy. Results for Different J_z could be extracted by multiplying E by J_xy.
            quantities = diagonalize_cluster_xxz(N=N, sectors_info=sectors_info, bonds=bonds, J_xy=1.0, J_z=J_z / J_xy)

            # printing diagonalization data to a file for save
            printing(quantities; name=diag_name, NTOP=topN, N=N, folder_path=diag_folder_path)

            # ---------------------------------------------------------

            # Here, we take the opportunity to read in the "lattice constants"
            # (multiplicities) for each topological cluster if we have reached
            # the end of the file.

            # skipping the subgraph information for now
            sc_Num = parse(Int64, readline(file))
            for s in 1:sc_Num
                readline(file)
            end
            readline(file)

            # Checking if we have reached the end of the file
            if occursin("#    Topo.    #    :    Multp. ", readline(file))
                # next(file)
                for i in 1:(topN+1)
                    multi[N, i] = parse(Int64, split(readline(file))[2])
                end
                EOF = true
                close(file)
                break
            else
                # If not yet done with clusters
                topN += 1
            end
        end
        println("finished diagonalizing clusters of order " * string(N))
    end
end


#diagonalize_all_clusters_xxz(5)
