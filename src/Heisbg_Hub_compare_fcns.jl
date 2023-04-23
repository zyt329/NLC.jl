using Plots
using JLD
using Dates
using Statistics
using Random, Distributions
using OffsetArrays: Origin
using OffsetArrays

"""
Does NLCE summation, return the raw summations.
"""
function NLCE(; Nmax, J, g, Temps, hs, multi)

    NT = length(Temps)
    Nh = length(hs)

    # O is going to hold NLCE's partial sums (sums involving clusters
    # in a particular order order) for our four properties. These are
    # Sn in Eq. (3) of https://arxiv.org/pdf/0706.3254.pdf
    # I think it needs to hold quantities for single sites at O[:,1,:]
    O = zeros(Float64, 6, Nmax + 1, NT, Nh)

    # This array is going to hold the contribution of clusters
    # from all orders
    weights = zeros(Float64, 6, Nmax + 1, 16000, NT, Nh)

    # Hard coding the contributions from the single site by hand
    # Be careful with the order of quantities in O and in julia code,
    # they might not be the same and caution is required.
    single_site_quants = single_site_quantities(Ts=Temps, hs=hs, g=g)
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
    #println("single site E(h=0.05):")
    #print(singleE[:,end])

    # Loop over NLCE orders
    for N = 2:Nmax
        # Opening the property binary files and loading the information
        # on to arrays
        floc = "../NLCE_Clusters/"
        fbase = floc * "NLCE_1_" * string(N)
        fname = fbase * ".txt"

        file = open(fname, "r")

        fname_quants = fbase * "_Heisbg" * "_quants.jld"
        quants_store = load(fname_quants)

        Estore = quants_store["Estore"]
        Mstore = quants_store["Mstore"]
        Esqstore = quants_store["Esqstore"]
        Msqstore = quants_store["Msqstore"]
        Nstore = quants_store["Nstore"]
        lnZstore = quants_store["lnZstore"]
        # print(Nstore[0][-1])
        # print("\n")

        # Skiping line with order number & toplogical number
        readline(file)

        # Tolopogy number
        c = parse(Int64, readline(file)) + 1 # +1 for julia indexing

        # Going through the cluster files again to read in the
        # subcluster information
        EOF = false
        while EOF == false

            line = split(readline(file))

            # Skiping the bond info
            for b = 1:parse(Int64, line[1])
                readline(file)
            end

            # Subgraph info
            sc_Num = parse(Int64, readline(file))
            for sg = 1:sc_Num
                line = split(readline(file))
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

            readline(file)
            if occursin("Topo", readline(file))
                # It's the end of the file
                EOF = true
                close(file)
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
    return O
end

"""
Do NLCE To the maximum order of "Nmax" and Wynn summation
for the given J,g at a series of temperatures "Temps", and
a series of magnetic fields "hs". (Only single magnetic field
because of the way the code's written.)

return: an np array of M values obtained by the Wynn summation
at the given "Temps" and "h"
"""
function M_Wynn_Euler(Nmax, J, g, Temps, hs, savepath, savename; nwynn=4, raw_sum_order=10)
    # Interaction J (in kelvin)
    # g factor

    # Holding "lattice constants" (multiplicities) for all clusters
    # in all orders. Probably good for up to order 9. Need to check
    # what is needed exactly for the size of the second dimension.
    multi = zeros(Int64, Nmax + 1, 16000)

    # load julia script
    # include("../Heisenberg.jl")
    # include("../measurements_Heisenberg.jl")

    name = "Heisbg_order=" * string(12) #string(Nmax) #always read from Nmax=12th order data 

    # Number of temperatures
    NT = length(Temps)

    # Number of magnetic field h
    Nh = length(hs)


    # Setting physical constants to 1.0
    #mu_B = 1.0  # joules per tesla
    #k_B = 1.0  # joules per K
    #N_A = 1.0  # Avogadro's number

    # Counter for isomorphically distinct clusters
    topN = 0
    # Loop over the NLCE order
    for N = 2:Nmax
        # Initializing arrays and openning files
        Estore = zeros(Float64, 16000, NT, Nh)
        Mstore = zeros(Float64, 16000, NT, Nh)
        Esqstore = zeros(Float64, 16000, NT, Nh)
        Msqstore = zeros(Float64, 16000, NT, Nh)
        Nstore = zeros(Float64, 16000, NT, Nh)
        lnZstore = zeros(Float64, 16000, NT, Nh)

        # Change 1 to 2 in the line below and in the next cell
        # to include n.n.n. bonds
        floc = "../NLCE_Clusters/"
        fbase = floc * "NLCE_1_" * string(N)
        fname = fbase * ".txt"
        file = open(fname, "r")
        fname_quants = fbase * "_Heisbg" * "_quants.jld"

        # Skips line with order number
        readline(file)

        # Going through the file for each order line by line
        # and reading in the cluster information
        topN = parse(Int64, readline(file))
        println("ORDER", N)

        EOF = false
        while EOF == false
            # NTOP starts from 1 for convenience of indexing in julia
            # topN starts from 0 as used in the cluster data from Ehsan
            NTOP = topN + 1
            line = split(readline(file))

            # Get the number of bonds from this line
            nB = parse(Int64, line[1])

            # Skip the lines with bond information
            for b = 1:nB
                line = split(readline(file))
                # next(file)
            end

            # Doing and storing quantieis:
            #
            # 1- Average energy, <H>
            # 2- Average magnetization, <M>
            # 3- Average energy squared, <H^2>, and
            # 4- Average magnetization squared, <M^2>
            # 5- Average total number of particle <N>

            # Thermal sums are done here **insert julia calculation here!!!
            # start = time.time()
            # Main.eval('cd("../")')
            quantities =
                reading_quantities(name=name, NTOP=topN, N=N, folder_path="../")

            quantities_avg = thermal_avg_loop(;
                Temps=Temps,
                hs=hs,
                J=J,
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
            # the end of the file.

            # skipping the subgraph information for now
            # print(line)
            sc_Num = parse(Int64, readline(file))
            for s = 1:sc_Num
                readline(file)
            end
            readline(file)

            # Checking if we have reached the end of the file
            if occursin("Topo", readline(file))
                # next(file)
                for i = 1:(topN+1)
                    multi[N, i] = parse(Int64, split(readline(file))[2])
                end
                EOF = true
                close(file)
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
            fname_quants,
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

    # In this part, we use the ED results, the subgraph information
    # and the lattice constants and construct the NLCE sums.

    # This array is going to hold the contribution of clusters
    # from all orders
    # weights = zeros(Float64, 6, Nmax + 1, 3500, NT, Nh)

    # Cluster counter; this is the same as topN from the previous part
    c = 0

    # sinclude("./Heisbg_Hub_compare_single_site.jl")



    O = NLCE(Nmax=Nmax, J=J, g=g, Temps=Temps, hs=hs, multi=multi)

    save(savepath * "O_" * savename * ".jld", "O", O)
    #println("O: ")
    #print(O[2,:,1,1])


    # ---------------------------------------------------------
    # In this part, we add partial sums to obtain the NLCE sums
    # upto different orders. We can do it without any tricks
    # (raw sums), or using numerical resummation tricks (here
    # we have implemented the Wynn algorithm, see Sec. 4 of
    # https://arxiv.org/pdf/1207.3366.pdf for more details).
    # The Euler resummation algorithm is also useful and worth
    # implementing.

    # Doing the raw sums for NLCE
    raw_O = zeros(Float64, 6, Nmax, NT, Nh)
    for N = 1:(Nmax)
        if N != 1
            raw_O[:, N, :, :] = raw_O[:, N-1, :, :] + O[:, N, :, :]
        else
            raw_O[:, N, :, :] = O[:, N, :, :]
        end
    end
    #println("raw sums: ")
    println("M raw sums:" * "$(raw_O[2, :, end, end])")


    # Doing the sums according to the Wynn algorithm for faster convergence
    # nwynn = 4
    #=
    epsln_O = zeros(Float64, 6, Nmax, NT, Nh, 16)
    for i = 1:NT
        for j = 1:Nh
            for N = 1:Nmax
                if N != 1
                    epsln_O[:, N, i, j, 1] = epsln_O[:, N-1, i, j, 1] + O[:, N, i, j]
                else
                    epsln_O[:, N, i, j, 1] = O[:, N, i, j]
                end
            end

            for k = 2:(2*nwynn+1)
                for N = 1:(Nmax-k+1)
                    for l = 1:6
                        delta = (epsln_O[l, N+1, i, j, k-1] - epsln_O[l, N, i, j, k-1])
                        # if abs(delta/epsln_O[j,N,i,k]) > 0.0:
                        if k == 2
                            epsln_O[l, N, i, j, k] = 1.0 / delta
                        else
                            epsln_O[l, N, i, j, k] =
                                (epsln_O[l, N+1, i, j, k-2] + 1.0 / delta)
                        end
                    end
                end
            end
        end
    end

    nwynn *= 2
    #println("wynn sums: ")
    #println(epsln_O[2, 5, 1, 1, 5])=#

    # try to write Wynn algorithm by shifting Array index to start from -1 by using package OffsetArrays
    # nwynn = 5
    epsln_O = Origin(1, 1, 1, 1, -1)(zeros(Float64, 6, Nmax, NT, Nh, 18))
    epsln_O[:, 1:Nmax, :, :, -1] = zeros(Float64, 6, Nmax, NT, Nh)
    epsln_O[:, 1:Nmax, :, :, 0] = raw_O[:, 1:Nmax, :, :]
    for k = 1:(2*nwynn)
        for n = 1:(Nmax-k)
            delta = epsln_O[:, n+1, :, :, k-1] - epsln_O[:, n, :, :, k-1]
            epsln_O[:, n, :, :, k] = (epsln_O[:, n+1, :, :, k-2] + 1 ./ delta)
        end
    end

    # Doing Euler summation
    # raw_sum_order = 10
    m_max = Nmax - raw_sum_order + 1
    Euler_O = zeros(Float64, 6, NT, Nh, m_max)
    # m=1 is the raw sum to the order of raw_sum_order
    Euler_O[:, :, :, 1] = raw_O[:, raw_sum_order, :, :]
    for m = 2:m_max
        Euler_O[:, :, :, m] = Euler_O[:, :, :, m-1]
        for n = 1:(m-1)
            Euler_O[:, :, :, m] +=
                (1 / 2)^(m - 1) * O[:, raw_sum_order+n, :, :] * binomial(m - 1, n)
        end
    end

    return [OffsetArrays.no_offset_view(epsln_O)[:, Nmax-2*nwynn, :, :, 2*nwynn+2], Euler_O[:, :, :, m_max], raw_O]
end
