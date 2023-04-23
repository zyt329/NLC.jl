using NLC
using Printf
using JLD

function run_resummation(; Nmax, Euler_raw_sum_order, nwynn, simulation_folder_full_path)

    #==========================================#
    #= READ IN SIMULATION INFO AND PARAMETERS =#
    #==========================================#

    # NLC result at all orders
    O = load(simulation_folder_full_path * "O.jld")["O"]

    # simulation info
    simulation_info = load(simulation_folder_full_path * "simulation_parameters.jld")

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

    #==================#
    #= DO RESUMMATION =#
    #==================#

    # store resummation info to a dict
    resummation_info = Dict("Nmax" => Nmax, "Euler_raw_sum_order" => Euler_raw_sum_order, "nwynn" => nwynn)

    # do the resummation
    resummation_data = resummation(; Nmax=Nmax, raw_sum_order=Euler_raw_sum_order, nwynn=nwynn, Ofilepath=simulation_folder_full_path, Ofilename="O")


    # load resummation data to arrays
    Wynn_data = resummation_data[1]
    Euler_data = resummation_data[2]
    raw_data = resummation_data[3][:, end, :, :]

    #=========================#
    #= Pick out M and S data =#
    #=========================#

    # M values
    M_Wynn = -Wynn_data[2, :, :]
    M_Euler = -Euler_data[2, :, :]
    M_raw = -raw_data[2, :, :]

    # extract S values
    S_Wynn = Wynn_data[6, :, :]
    S_Euler = Euler_data[6, :, :]
    S_raw = raw_data[6, :, :]
    for (T_ind, T) in enumerate(Temps)
        for (h_ind, h) in enumerate(hs)
            S_Wynn[T_ind, h_ind] +=
                1 / T * Wynn_data[1, T_ind, h_ind]
            S_Euler[T_ind, h_ind] +=
                1 / T * Euler_data[1, T_ind, h_ind]
            S_raw[T_ind, h_ind] += 1 / T * raw_data[1, T_ind, h_ind]
        end
    end

    #=================================#
    #= SAVE RESUMMATION DATA TO FILE =#
    #=================================#

    open(resummation_folder_path * "/resummation_data.csv", "w") do io

        # write file header
        write(
            io,
            "T h M_Wynn S_Wynn M_Euler S_Euler M_raw S_raw \n",
        )

        # write data to file
        for (h_ind, h) in enumerate(hs)
            for (T_ind, T) in enumerate(Temps)
                write(io, string(T) * " ")
                write(io, string(h) * " ")
                write(io, string(M_Wynn[T_ind, h_ind]) * " ")
                write(io, string(S_Wynn[T_ind, h_ind]) * " ")
                write(io, string(M_Euler[T_ind, h_ind]) * " ")
                write(io, string(S_Euler[T_ind, h_ind]) * " ")
                write(io, string(M_raw[T_ind, h_ind]) * " ")
                write(io, string(S_raw[T_ind, h_ind]) * " ")
                write(io, "\n")
            end
        end

    end # done writing

    # save to .jld file
    save(
        resummation_folder_path * "/resummation_data.jld",
        "Temps",
        Temps,
        "hs",
        hs,
        "M_Wynn",
        M_Wynn,
        "S_Wynn",
        S_Wynn,
        "M_Euler",
        M_Euler,
        "S_Euler",
        S_Euler,
        "M_raw",
        M_raw,
        "S_raw",
        S_raw
    )

    print("Resummation completed")
end


run_resummation(Nmax=9, Euler_raw_sum_order=4, nwynn=2, simulation_folder_full_path="/nfs/home/zyt329/NLC/examples/simulation_J_z[J_xy1.0000-1/")
