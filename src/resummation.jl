#using Plots
using JLD
using Dates
#using Statistics
#using Random, Distributions
using OffsetArrays:Origin
using OffsetArrays


Ofilepath = "E:/UC Davis/Research/triangular_Hubbard/nlceforhubbardontriangularlattice/Heisbg_Hub_compare/results/"
Ofilename = "O_Heisbg_few_T_sweep_h_result_20_Nov_2022_21_10"

"""
This script prints the resummation data to a .txt file
"""
function print_data(;Nmax=12, raw_sum_order1=4, raw_sum_order2=4, nwynn=5)
    cd("E:/UC Davis/Research/triangular_Hubbard/nlceforhubbardontriangularlattice/Heisbg_Hub_compare/")
    # write result to a file
    my_time = Dates.now()
    savepath = Ofilepath
    savename = "Heisbg_few_T_sweep_h_result_$(Dates.format(my_time, "dd_u_yyyy_HH_MM"))_nwynn$(nwynn)_Euler$(raw_sum_order1)"
    open(savepath * savename * ".txt", "w") do io
        # Nmax from input argument
        # Choosing a theoretical logarithmic temperature grid
        # from 0.7K to 350K
        Temps = [0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9]
        # Temps = vcat(range(0.05, 10, length = NT))
        hs = vcat(range(0.0, 1, length = 100))
        # hs = exp10.(range(-1, 1, length=100))

        # Number of temperatures
        NT = length(Temps)
        # Number of fields h
        Nh = length(hs)

        # calculate NLC and Wynn sums
        resum_data_9order = resummation(;Nmax=Nmax, raw_sum_order=raw_sum_order1, nwynn = nwynn, Ofilepath=Ofilepath, Ofilename=Ofilename)
        resum_data_8order = resummation(;Nmax=Nmax-1, raw_sum_order=raw_sum_order2, nwynn = nwynn, Ofilepath=Ofilepath, Ofilename=Ofilename)
        # resum_data_8order = resum_data_9order #redundant, place holder (懒得改整篇代码了)
        Wynn_data_9order = resum_data_9order[1]
        Wynn_data_8order = resum_data_8order[1]
        Euler_data_9order = resum_data_9order[2]
        Euler_data_8order = resum_data_8order[2]
        raw_data_9order = resum_data_9order[3][:, end, :, :]
        raw_data_8order = resum_data_8order[3][:, end, :, :]

        # extract M values
        M_Wynn_data_9order = -Wynn_data_9order[2, :, :]
        M_Wynn_data_8order = -Wynn_data_8order[2, :, :]
        M_Euler_data_9order = -Euler_data_9order[2, :, :]
        M_Euler_data_8order = -Euler_data_8order[2, :, :]
        M_raw_data_9order = -raw_data_9order[2, :, :]
        M_raw_data_8order = -raw_data_8order[2, :, :]
        # extract S values
        S_Wynn_data_9order = Wynn_data_9order[6, :, :]
        S_Wynn_data_8order = Wynn_data_8order[6, :, :]
        S_Euler_data_9order = Euler_data_9order[6, :, :]
        S_Euler_data_8order = Euler_data_8order[6, :, :]
        S_raw_data_9order = raw_data_9order[6, :, :]
        S_raw_data_8order = raw_data_8order[6, :, :]
        for (T_ind, T) in enumerate(Temps)
            for (h_ind, h) in enumerate(hs)
                S_Wynn_data_9order[T_ind, h_ind] +=
                    1 / T * Wynn_data_9order[1, T_ind, h_ind]
                S_Wynn_data_8order[T_ind, h_ind] +=
                    1 / T * Wynn_data_8order[1, T_ind, h_ind]
                S_Euler_data_9order[T_ind, h_ind] +=
                    1 / T * Euler_data_9order[1, T_ind, h_ind]
                S_Euler_data_8order[T_ind, h_ind] +=
                    1 / T * Euler_data_8order[1, T_ind, h_ind]
                S_raw_data_9order[T_ind, h_ind] += 1 / T * raw_data_9order[1, T_ind, h_ind]
                S_raw_data_8order[T_ind, h_ind] += 1 / T * raw_data_8order[1, T_ind, h_ind]
            end
        end

        # calculate distance btw experiment and theory for M for this set of parameter.
        write(
            io,
            "T h M$(Nmax) M$(Nmax-1) S$(Nmax) S$(Nmax-1) M_Euler$(Nmax) M_Euler$(Nmax-1) S_Euler$(Nmax) S_Euler$(Nmax-1) M_raw$(Nmax) M_raw$(Nmax-1) S_raw$(Nmax) S_raw$(Nmax-1)\n",
        )
        for (h_ind, h) in enumerate(hs)
            for (T_ind, T) in enumerate(Temps)
                write(io, string(T) * " ")
                write(io, string(h) * " ")
                write(io, string(M_Wynn_data_9order[T_ind, h_ind]) * " ")
                write(io, string(M_Wynn_data_8order[T_ind, h_ind]) * " ")
                write(io, string(S_Wynn_data_9order[T_ind, h_ind]) * " ")
                write(io, string(S_Wynn_data_8order[T_ind, h_ind]) * " ")
                write(io, string(M_Euler_data_9order[T_ind, h_ind]) * " ")
                write(io, string(M_Euler_data_8order[T_ind, h_ind]) * " ")
                write(io, string(S_Euler_data_9order[T_ind, h_ind]) * " ")
                write(io, string(S_Euler_data_8order[T_ind, h_ind]) * " ")
                write(io, string(M_raw_data_9order[T_ind, h_ind]) * " ")
                write(io, string(M_raw_data_8order[T_ind, h_ind]) * " ")
                write(io, string(S_raw_data_9order[T_ind, h_ind]) * " ")
                write(io, string(S_raw_data_8order[T_ind, h_ind]) * " \n")
            end
        end

        # save to .jld files
        save(
            savepath * savename * ".jld",
            "M$(Nmax)",
            M_Wynn_data_9order,
            "M$(Nmax-1)",
            M_Wynn_data_8order,
            "Temps",
            Temps,
            "hs",
            hs,
            "S$(Nmax)",
            S_Wynn_data_9order,
            "S$(Nmax-1)",
            S_Wynn_data_8order,
            "M_Euler$(Nmax)",
            M_Euler_data_9order,
            "M_Euler$(Nmax-1)",
            M_Euler_data_8order,
            "S_Euler$(Nmax)",
            S_Euler_data_9order,
            "S_Euler$(Nmax-1)",
            S_Euler_data_8order,
            "M_raw$(Nmax)",
            M_raw_data_9order,
            "M_raw$(Nmax-1)",
            M_raw_data_8order,
            "S_raw$(Nmax)",
            S_raw_data_9order,
            "S_raw$(Nmax-1)",
            S_raw_data_8order,
        )

        print("finished calculation.")
    end
    #return the full path to the result file.
    return savepath * savename * ".jld"
end



function resummation(;Nmax, raw_sum_order, nwynn, Ofilepath, Ofilename)

    O = load(Ofilepath * Ofilename * ".jld")["O"]
    O_size = size(O)
    Nmax_raw = O_size[2]
    NT = O_size[3]
    Nh = O_size[4]
    @assert(Nmax <= Nmax_raw, "maximum order must be smaller than that of the raw data: " * string(Nmax_raw))
    @assert(raw_sum_order <= Nmax, "Euler sum order must be smaller than Nmax: " * string(Nmax))
    # ---------------------------------------------------------
    # O = zeros(Float64, 6, Nmax + 1, NT, Nh)
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
    #print(raw_O[2,:,1,1])


    # try to write Wynn algorithm by shifting Array index to start from -1 by using package OffsetArrays
    # nwynn = 5
    epsln_O = Origin(1,1,1,1,-1)(zeros(Float64, 6, Nmax, NT, Nh, 18))
    epsln_O[:,1:Nmax,:,:,-1] = zeros(Float64, 6, Nmax, NT, Nh)
    epsln_O[:,1:Nmax,:,:,0] = raw_O[:, 1:Nmax, :, :]
    for k = 1:(2*nwynn)
        for n = 1:(Nmax-k)
            delta = epsln_O[:, n+1, :, :, k-1] - epsln_O[:, n, :, :, k-1]
            epsln_O[:, n, :, :, k] = (epsln_O[:, n+1, :, :, k-2] + 1 ./ delta )
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

    # return [epsln_O[:, Nmax-nwynn, :, :, nwynn+1], Euler_O[:, :, :, m_max], raw_O]
    
    return [OffsetArrays.no_offset_view(epsln_O)[:, Nmax-2*nwynn, :, :, 2*nwynn+2], Euler_O[:, :, :, m_max], raw_O]
end

print_data()