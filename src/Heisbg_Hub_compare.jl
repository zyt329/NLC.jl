using DelimitedFiles
using Dates
using JLD

cd("E:/UC Davis/Research/triangular_Hubbard/nlceforhubbardontriangularlattice/Heisbg_Hub_compare/")
include("./Heisbg_Hub_compare_fcns.jl")
include("./measurements_Heisbg_Hub_compare.jl")
include("./Heisbg_Hub_compare_single_site.jl")

"""
This script loops over selected (Temp, h) values. Calculate M and S for comparasion with Hubbard model data from Owen.
"""
function driver(Nmax)
    
    # write result to a file
    my_time = Dates.now()
    savepath = "./results/M_vs_T/"
    savename = "Heisbg_few_h1_3_sweep_T_result_$(Dates.format(my_time, "dd_u_yyyy_HH_MM"))"
    open(savepath * savename * ".txt", "w") do io
        # Choosing a theoretical logarithmic temperature grid
        # from 0.7K to 350K
        Temps = vcat(range(0.0, 1, length = 100))
        # Temps = vcat(range(0.05, 10, length = NT))
        hs = vcat(range(1.0, 3.0, length = 11))
        # hs = exp10.(range(-1, 1, length=100))

        # Number of temperatures
        NT = length(Temps)
        # Number of fields h
        Nh = length(hs)

        # Setting all physical constants to 1.0
        J1 = 1.0
        g1 = 1.0
        mu_B = 1.0  # joules per tesla
        k_B = 1.0  # joules per K
        N_A = 1.0  # Avogadro's number

        # calculate NLC and Wynn sums
        resum_data_9order = M_Wynn_Euler(Nmax, J1, g1, Temps, hs, savepath, savename)
        #resum_data_8order = M_Wynn_Euler(Nmax, J1, g1, Temps, hs)
        resum_data_8order = resum_data_9order #redundant, place holder (懒得改整篇代码了)
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
    return nothing
end

@time driver(12)
