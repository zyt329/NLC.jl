using NLC
using Printf
using JLD2

function run_resummation(; Nmax_data, simulation_folder_full_path)

    #==========================================#
    #= READ IN SIMULATION INFO AND PARAMETERS =#
    #==========================================#

    # NLC result at all orders
    O = load(simulation_folder_full_path * "O.jld2")["O"]

    # simulation info
    simulation_info = load(simulation_folder_full_path * "simulation_parameters.jld2")

    Temps = simulation_info["Temps"]
    NT = length(Temps)

    hs = simulation_info["hs"]
    Nh = length(hs)

    J_xy = simulation_info["J_xy"]

    J_z = simulation_info["J_z"]

    k_B = simulation_info["k_B"]

    N_A = simulation_info["N_A"]

    mu_B = simulation_info["mu_B"]

    g = simulation_info["g"]

    #========================================#
    #= MAKE FOLDER TO HOLD RESUMMATION DATA =#
    #========================================#

    resummation_folder_path = make_indexed_folder(folder_prefix="resummation", folder_path=simulation_folder_full_path)

    # doing resummation for all Nmax and allowed Euler_raw_sum_order, nwynn
    for Nmax in 4:Nmax_data

        println("resumming Nmax=$Nmax")

        # all orders of Euler resummation
        for Euler_raw_sum_order in 1:(Nmax-1)

            # do the resummation
            resummation_data = resummation(; Nmax=Nmax, raw_sum_order=Euler_raw_sum_order, nwynn=1, Ofilepath=simulation_folder_full_path, Ofilename="O")

            Euler_data = resummation_data[2]

            M_Euler = -Euler_data[2, :, :]
            E_Euler = Euler_data[1, :, :]
            S_Euler = Euler_data[6, :, :]
            # extract S values
            for (T_ind, T) in enumerate(Temps)
                for (h_ind, h) in enumerate(hs)
                    S_Euler[T_ind, h_ind] +=
                        1 / T * Euler_data[1, T_ind, h_ind]
                end
            end
            # extract Chi values
            Chi_Euler = Euler_data[4, :, :]
            for (T_ind, T) in enumerate(Temps)
                for (h_ind, h) in enumerate(hs)
                    Chi_Euler[T_ind, h_ind] =
                        1 / T * Chi_Euler[T_ind, h_ind]
                end
            end
            # extract C values
            C_Euler = Euler_data[3, :, :]
            for (T_ind, T) in enumerate(Temps)
                for (h_ind, h) in enumerate(hs)
                    C_Euler[T_ind, h_ind] =
                        1 / T^2 * C_Euler[T_ind, h_ind]
                end
            end

            # write to csv
            Euler_data_name = @sprintf "/resummation_data_Nmax%d_Euler%d" Nmax Euler_raw_sum_order
            open(resummation_folder_path * Euler_data_name * ".csv", "w") do io

                # write file header
                write(
                    io,
                    "T h M S Chi C E \n",
                )

                # write data to file
                for (h_ind, h) in enumerate(hs)
                    for (T_ind, T) in enumerate(Temps)
                        write(io, string(T) * " ")
                        write(io, string(h) * " ")
                        write(io, string(M_Euler[T_ind, h_ind]) * " ")
                        write(io, string(S_Euler[T_ind, h_ind]) * " ")
                        write(io, string(Chi_Euler[T_ind, h_ind]) * " ")
                        write(io, string(C_Euler[T_ind, h_ind]) * " ")
                        write(io, string(E_Euler[T_ind, h_ind]) * " ")
                        write(io, "\n")
                    end
                end
            end # done writing

            # save to .jld file
            save(
                resummation_folder_path * Euler_data_name * ".jld2",
                "Temps",
                Temps,
                "hs",
                hs,
                "M",
                M_Euler,
                "S",
                S_Euler,
                "Chi",
                Chi_Euler,
                "C",
                C_Euler,
                "E",
                E_Euler,
                "Euler_raw_sum_order",
                Euler_raw_sum_order,
                "Nmax",
                Nmax
            )


        end

        # all orders of Wynn resummation
        for nwynn in 1:(Int(floor((Nmax - 1) / 2)))

            # do the resummation
            resummation_data = resummation(; Nmax=Nmax, raw_sum_order=(Nmax - 1), nwynn=nwynn, Ofilepath=simulation_folder_full_path, Ofilename="O")

            Wynn_data = resummation_data[1]

            M_Wynn = -Wynn_data[2, :, :]
            E_Wynn = Wynn_data[1, :, :]
            S_Wynn = Wynn_data[6, :, :]
            # extract S values
            for (T_ind, T) in enumerate(Temps)
                for (h_ind, h) in enumerate(hs)
                    S_Wynn[T_ind, h_ind] +=
                        1 / T * Wynn_data[1, T_ind, h_ind]
                end
            end
            # extract Chi values
            Chi_Wynn = Wynn_data[4, :, :]
            for (T_ind, T) in enumerate(Temps)
                for (h_ind, h) in enumerate(hs)
                    Chi_Wynn[T_ind, h_ind] =
                        1 / T * Chi_Wynn[T_ind, h_ind]
                end
            end
            # extract C values
            C_Wynn = Wynn_data[3, :, :]
            for (T_ind, T) in enumerate(Temps)
                for (h_ind, h) in enumerate(hs)
                    C_Wynn[T_ind, h_ind] =
                        1 / T^2 * C_Wynn[T_ind, h_ind]
                end
            end

            # write to csv
            Wynn_data_name = @sprintf "/resummation_data_Nmax%d_Wynn%d" Nmax nwynn
            open(resummation_folder_path * Wynn_data_name * ".csv", "w") do io

                # write file header
                write(
                    io,
                    "T h M S Chi C E \n",
                )

                # write data to file
                for (h_ind, h) in enumerate(hs)
                    for (T_ind, T) in enumerate(Temps)
                        write(io, string(T) * " ")
                        write(io, string(h) * " ")
                        write(io, string(M_Wynn[T_ind, h_ind]) * " ")
                        write(io, string(S_Wynn[T_ind, h_ind]) * " ")
                        write(io, string(Chi_Wynn[T_ind, h_ind]) * " ")
                        write(io, string(C_Wynn[T_ind, h_ind]) * " ")
                        write(io, string(E_Wynn[T_ind, h_ind]) * " ")
                        write(io, "\n")
                    end
                end
            end # done writing

            # save to .jld file
            save(
                resummation_folder_path * Wynn_data_name * ".jld2",
                "Temps",
                Temps,
                "hs",
                hs,
                "M",
                M_Wynn,
                "S",
                S_Wynn,
                "Chi",
                Chi_Wynn,
                "C",
                C_Wynn,
                "E",
                E_Wynn,
                "nwynn",
                nwynn,
                "Nmax",
                Nmax
            )
        end

        # ========================== #
        # = all orders of raw data = #
        # ========================== #

        # do the resummation
        resummation_data = resummation(; Nmax=Nmax, raw_sum_order=(Nmax - 1), nwynn=1, Ofilepath=simulation_folder_full_path, Ofilename="O")

        raw_data = resummation_data[3][:, end, :, :]

        M_raw = -raw_data[2, :, :]
        E_raw = raw_data[1, :, :]
        S_raw = raw_data[6, :, :]
        # extract S values
        for (T_ind, T) in enumerate(Temps)
            for (h_ind, h) in enumerate(hs)
                S_raw[T_ind, h_ind] +=
                    1 / T * raw_data[1, T_ind, h_ind]
            end
        end
        # extract Chi values
        Chi_raw = raw_data[4, :, :]
        for (T_ind, T) in enumerate(Temps)
            for (h_ind, h) in enumerate(hs)
                Chi_raw[T_ind, h_ind] =
                    1 / T * Chi_raw[T_ind, h_ind]
            end
        end
        # extract C values
        C_raw = raw_data[3, :, :]
        for (T_ind, T) in enumerate(Temps)
            for (h_ind, h) in enumerate(hs)
                C_raw[T_ind, h_ind] =
                    1 / T^2 * C_raw[T_ind, h_ind]
            end
        end

        # write to csv
        raw_data_name = @sprintf "/raw_data_Nmax%d" Nmax
        open(resummation_folder_path * raw_data_name * ".csv", "w") do io

            # write file header
            write(
                io,
                "T h M S Chi C E \n",
            )

            # write data to file
            for (h_ind, h) in enumerate(hs)
                for (T_ind, T) in enumerate(Temps)
                    write(io, string(T) * " ")
                    write(io, string(h) * " ")
                    write(io, string(M_raw[T_ind, h_ind]) * " ")
                    write(io, string(S_raw[T_ind, h_ind]) * " ")
                    write(io, string(Chi_raw[T_ind, h_ind]) * " ")
                    write(io, string(C_raw[T_ind, h_ind]) * " ")
                    write(io, string(E_raw[T_ind, h_ind]) * " ")
                    write(io, "\n")
                end
            end
        end # done writing

        # save to .jld file
        save(
            resummation_folder_path * raw_data_name * ".jld2",
            "Temps",
            Temps,
            "hs",
            hs,
            "M",
            M_raw,
            "S",
            S_raw,
            "Chi",
            Chi_raw,
            "C",
            C_raw,
            "E",
            E_raw,
            "Nmax",
            Nmax
        )


    end

    println("Resummation completed")
end


run_resummation(Nmax_data=14, simulation_folder_full_path="/nfs/home/zyt329/Research/xxz/runs/hashtag_simulation_J_z[J_xy32.0000-2/")
