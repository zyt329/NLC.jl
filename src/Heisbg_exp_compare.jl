using DelimitedFiles
using Dates
using JLD

cd("E:/UC Davis/Research/triangular_Hubbard/nlceforhubbardontriangularlattice/Heisbg_Hub_compare/")
include("./consts.jl")
include("./Heisbg_Hub_compare_fcns.jl")
include("./measurements_Heisbg_Hub_compare.jl")
include("./Heisbg_Hub_compare_single_site.jl")

"""
This script loops over use the chosen (J1,J2,g1,g2) and compare M with experiment.
"""
function driver(;Nmax=12, J1=11.75,g1=2.075,J2=3.25,g2=2.15)
    # write result to a file
    my_time = Dates.now()
    savepath = "./results/Heisbg_exp_compare/"
    savename = "Heisbg_exp_compare_$(Dates.format(my_time, "dd_u_yyyy_HH_MM"))"
    open(savepath * savename * ".txt", "w") do io
        # Choosing a theoretical logarithmic temperature grid
        # from 0.7K to 350K
        Temps = exp10.(vcat(range(log10(0.4*J1), log10(300), length = 50)))
        # Temps = vcat(range(0.05, 10, length = NT))
        hs = [1, 3, 5, 6, 7]
        # hs = exp10.(range(-1, 1, length=100))

        # Number of temperatures
        NT = length(Temps)
        # Number of fields h
        Nh = length(hs)

        # calculate NLC and Wynn sums
        resum_data_1 = M_Wynn_Euler(Nmax, J1, g1, Temps, hs, savepath, savename) .* N_A
        resum_data_2 = M_Wynn_Euler(Nmax, J2, g2, Temps, hs, savepath, savename) .* N_A
        Wynn_data_1 = resum_data_1[1]
        Wynn_data_2 = resum_data_2[1]
        Wynn_data = (Wynn_data_1 + Wynn_data_2) ./ 2
        Euler_data_1 = resum_data_1[2]
        Euler_data_2 = resum_data_2[2]
        Euler_data = (Euler_data_1 + Euler_data_2) ./ 2
        raw_data_1 = resum_data_1[3][:, end, :, :]
        raw_data_2 = resum_data_2[3][:, end, :, :]
        raw_data = (raw_data_1 +  raw_data_2) ./ 2

        # extract M values
        M_Wynn_data = -Wynn_data[2, :, :]
        M_Euler_data = -Euler_data[2, :, :]
        M_raw_data = -raw_data[2, :, :]
        # extract S values
        S_Wynn_data = Wynn_data[6, :, :]
        S_Euler_data = Euler_data[6, :, :]
        S_raw_data = raw_data[6, :, :]
        for (T_ind, T) in enumerate(Temps)
            for (h_ind, h) in enumerate(hs)
                S_Wynn_data[T_ind, h_ind] +=
                    1 / T * Wynn_data[1, T_ind, h_ind]
                S_Wynn_data[T_ind, h_ind] +=
                    1 / T * Wynn_data[1, T_ind, h_ind]
                S_Euler_data[T_ind, h_ind] +=
                    1 / T * Euler_data_1[1, T_ind, h_ind]
                S_Euler_data[T_ind, h_ind] +=
                    1 / T * Euler_data[1, T_ind, h_ind]
                S_raw_data[T_ind, h_ind] += 1 / T * raw_data[1, T_ind, h_ind]
                S_raw_data[T_ind, h_ind] += 1 / T * raw_data[1, T_ind, h_ind]
            end
        end

        # calculate distance btw experiment and theory for M for this set of parameter.
        write(
            io,
            "T h M$(Nmax) S$(Nmax) M_Euler$(Nmax) S_Euler$(Nmax) M_raw$(Nmax) S_raw$(Nmax) \n",
        )
        for (h_ind, h) in enumerate(hs)
            for (T_ind, T) in enumerate(Temps)
                write(io, string(T) * " ")
                write(io, string(h) * " ")
                write(io, string(M_Wynn_data[T_ind, h_ind]) * " ")
                write(io, string(S_Wynn_data[T_ind, h_ind]) * " ")
                write(io, string(M_Euler_data[T_ind, h_ind]) * " ")
                write(io, string(S_Euler_data[T_ind, h_ind]) * " ")
                write(io, string(M_raw_data[T_ind, h_ind]) * " ")
                write(io, string(S_raw_data[T_ind, h_ind]) * " \n")
            end
        end



        # save to .jld files
        save(
            savepath * savename * ".jld",
            "M$(Nmax)",
            M_Wynn_data,
            "Temps",
            Temps,
            "hs",
            hs,
            "S$(Nmax)",
            S_Wynn_data,
            "M_Euler$(Nmax)",
            M_Euler_data,
            "S_Euler$(Nmax)",
            S_Euler_data,
            "M_raw$(Nmax)",
            M_raw_data,
            "S_raw$(Nmax)",
            S_raw_data,
        )

        print("finished calculation.")
    end
    return nothing
end

@time driver()
