#= functions to
    1. Diagonalize hamiltonian for the sectors of each cluster
    2. For each eigenstate, calculate expectation values of desired quantities
    3. Store the eigenvalues, <quantities> of each cluster in a file
    4. For each parameter(T,μ,h...), calculate Thermal average of desired quantities and feed back to the NLCE Python script.
=#
"""
    Function to produce eigenvalues(energies) of a cluster with J=1 and corresponding expectation values of quantities of interest. (eigen energies of other J values can be obtained by scaling the energy by En(J=J)=En(J=1)*J)
    Input:
        N::Float64 : number of sites in the cluster.
        sectors_info::Dict{Symbol,Any} : information of symmetry sectors of the cluster.
        bonds::Array{Array{Int,1},1} : bond information of the cluster.
    Output:
        [E, Quantity1, Quantity2...]. E, Quantity1, Quantity2... are all arrays of length 2^N
"""
function diagonalize_cluster_xxz(; N::Int64, sectors_info::Dict{Symbol,Any}, m0_sectors_info=nothing, bonds::Array{Array{Int,1},1}, J_xy::Float64=1.0, J_z::Float64)
    #P::Float64 = 0
    E = Float64[]
    Esq = Float64[]
    M = Float64[]
    Msq = Float64[]
    N_tot = Float64[] # place holder - to be deleted

    # loop over all symmetry sectors of number of up spins m           
    # only loop over half of the sectors
    # other sectors related by flipping all spins (ΠᵢSˣᵢ)
    for m = 0:Int(ceil(N / 2) - 1)

        # construct hamiltonian for the sector and diagonalize
        H_m = H_sector(; J_xy=J_xy, J_z=J_z, N=N, m=m, sectors_info=sectors_info, bonds=bonds)

        # diagonalize the sector
        (Es, states) = eigen(Hermitian(Matrix(H_m)))

        # calculate and record quantities for each eigenstate
        for (ind, En) in enumerate(Es)
            push!(E, En)
            push!(Esq, En^2)
            push!(M, (2m - N) * 1 / 2)
            push!(Msq, ((2m - N) * 1 / 2)^2)
            push!(N_tot, (N)) # place holder - to be deleted
        end

        # record quantities for Sᶻ=-m sector
        for (ind, En) in enumerate(Es)
            push!(E, En)
            push!(Esq, En^2)
            push!(M, (N - 2m) * 1 / 2)
            push!(Msq, ((N - 2m) * 1 / 2)^2)
            push!(N_tot, (N)) # place holder - to be deleted
        end
    end

    # treat Sᶻ=0 sector separately (when N is even)
    if iseven(N)

        # construct even and odd sectors' Hamiltonian
        H_m0_even, H_m0_odd = H_m0_sector(J_xy=J_xy, J_z=J_z, N=N, m0_sectors_info=m0_sectors_info, bonds=bonds)

        # diagonalize the even sector
        (Es, states) = eigen(Hermitian(Matrix(H_m0_even)))

        # calculate and record quantities for each eigenstate
        for (ind, En) in enumerate(Es)
            push!(E, En)
            push!(Esq, En^2)
            push!(M, 0)
            push!(Msq, 0)
            push!(N_tot, (N)) # place holder - to be deleted
        end

        # diagonalize the odd sector
        (Es, states) = eigen(Hermitian(Matrix(H_m0_odd)))

        # calculate and record quantities for each eigenstate
        for (ind, En) in enumerate(Es)
            push!(E, En)
            push!(Esq, En^2)
            push!(M, 0)
            push!(Msq, 0)
            push!(N_tot, (N)) # place holder - to be deleted
        end
    end

    return [E, Esq, M, Msq, N_tot]
end

"""
    do thermal average of quantities read from the file created by E_Quants() and printing()
    Input:
        J::Real : coupling in unit of Kelvin
        h::Real : magnetic field in unit of Tesla
        g::Real : g factor, pure number
    return:
        An array of thermal average of quantities at temperature of T
"""
function thermal_avg(; T::Real, J::Real, eig_vals, h::Real=0.0, g::Real, Ising=false)

    # passing in quantities
    if !Ising
        E::Array{Float64,1} = read(eig_vals["E"])
        M::Array{Float64,1} = read(eig_vals["M"])
        N_tot::Array{Float64,1} = ones(size(E)) # place holder
    else
        E = eig_vals["E"]
        M = eig_vals["M"]
        N_tot = ones(size(E)) # place holder
    end

    # calculate thermal average
    β = 1 / T

    # using matrix multiplication for speed
    E_with_h::Array{Float64,1} = (J * k_B) .* E .+ (h * g * mu_B) .* M
    P::Array{Float64,1} = exp.((-β / k_B) .* E_with_h)
    Z::Float64 = sum(P)
    E_avg::Float64 = sum(E_with_h .* P) / Z
    M_avg::Float64 = g * mu_B * sum(M .* P) / Z
    Esq_avg::Float64 = sum(E_with_h .^ 2 .* P) / Z
    Msq_avg::Float64 = (g * mu_B)^2 * sum(M .^ 2 .* P) / Z
    N_tot_avg::Float64 = sum(N_tot .* P) / Z


    """
    Z::Float64 = 0
    E_avg::Float64 = 0
    Esq_avg::Float64 = 0
    M_avg::Float64 = 0
    Msq_avg::Float64 = 0
    N_tot_avg::Float64 = 0
    for (n, En) in enumerate(E)
        P = exp(-β * (En * J + h * g * mu_B * M[n] / k_B))# shift energy?
        Z += P
        #E_avg += En*J * P  #doesn't include magnetic energy
        E_avg += (En * J * k_B + h * g * mu_B * M[n]) * P #include magnetic energy
        #Esq_avg += Esq[n]*J^2 * P # doesn't include magnetic energy
        Esq_avg += (En * J * k_B + h * g * mu_B * M[n])^2 * P
        M_avg += M[n] * g * mu_B * P
        Msq_avg += Msq[n] * (g * mu_B)^2 * P
        N_tot_avg += N_tot[n] * P
    end
    """
    return Dict("Z" => Z, "E" => E_avg, "Esq" => Esq_avg, "M" => M_avg, "Msq" => Msq_avg, "N" => N_tot_avg)
end


quant_names = ['E' "Esq" 'M' "Msq" "N_tot"]

"""
    Read in quantities(eigen energies and corresponding average quantities of eigenstates), a range of temperatures at which we want to calculate the thermal average. Out put thermal average of quantities at the read-in temperatures.

    Method for when h is passed in.
"""
function thermal_avg_hT_loop(; Temps, J::Real, eig_vals, hs::Vector{Type}=[0.0], g::Real=2.1, Ising=false) where {Type<:Real}
    # read in eigen values for the cluster (using HDF5)
    # eig_vals = h5open(diag_file_path, "r")

    NT = length(Temps)
    Nh = length(hs)
    Zs = Array{Float64}(undef, NT, Nh)
    E_avgs = Array{Float64}(undef, NT, Nh)
    Esq_avgs = Array{Float64}(undef, NT, Nh)
    M_avgs = Array{Float64}(undef, NT, Nh)
    Msq_avgs = Array{Float64}(undef, NT, Nh)
    N_tot_avgs = Array{Float64}(undef, NT, Nh)
    for (h_ind, h) in enumerate(hs)
        for (T_ind, T) in enumerate(Temps)
            avgs_T = thermal_avg(T=T, J=J, eig_vals=eig_vals, h=h, g=g, Ising=Ising)
            Zs[T_ind, h_ind] = avgs_T["Z"]
            E_avgs[T_ind, h_ind] = avgs_T["E"]
            Esq_avgs[T_ind, h_ind] = avgs_T["Esq"]
            M_avgs[T_ind, h_ind] = avgs_T["M"]
            Msq_avgs[T_ind, h_ind] = avgs_T["Msq"]
            N_tot_avgs[T_ind, h_ind] = avgs_T["N"]
        end
    end

    return Dict("Z" => Zs, "E" => E_avgs, "Esq" => Esq_avgs, "M" => M_avgs, "Msq" => Msq_avgs, "N" => N_tot_avgs)
end

"""
function thermal_avg_hT_loop(; Temps, J::Real, diag_file_path::String, hs::Vector{Type}=[0.0], g::Real=2.1, t_load, t_sum) where {Type<:Real}
    # read in eigen values for the cluster
    dt_load = @elapsed begin
        eig_vals = jldopen(diag_file_path, "r")
    end
    t_load[1] += dt_load

    NT = length(Temps)
    Nh = length(hs)

    Zs = Array{Float64}(undef, Nh, NT)
    E_avgs = Array{Float64}(undef, Nh, NT)
    Esq_avgs = Array{Float64}(undef, Nh, NT)
    M_avgs = Array{Float64}(undef, Nh, NT)
    Msq_avgs = Array{Float64}(undef, Nh, NT)
    N_tot_avgs = Array{Float64}(undef, Nh, NT)
    
    # expand every quantities to 3-d Matrix, Q[eig_val_index, h_index, T_index]
    E_with_h::Array{Float64,2} = ((J * k_B) * E)  + (h * g * mu_B) * M


    #=
    dt_sum = @elapsed begin
        # using matrix multiplication for speed
        E_with_h::Array{Float64,1} = (J * k_B) * E + (h * g * mu_B) * M
        P::Array{Float64,1} = exp.((-β / k_B) * E_with_h)
        Z::Float64 = sum(P)
        E_avg::Float64 = sum(E_with_h .* P) / Z
        M_avg::Float64 = g * mu_B * sum(M .* P) / Z
        Esq_avg::Float64 = sum(E_with_h .^ 2 .* P) / Z
        Msq_avg::Float64 = (g * mu_B)^2 * sum(M .^ 2 .* P) / Z
        N_tot_avg::Float64 = sum(N_tot .* P) / Z
    end


    =#



    avgs_T = thermal_avg(T=T, J=J, eig_vals=eig_vals, h=h, g=g, t_load=t_load, t_sum=t_sum)
    Zs[T_ind, h_ind] = avgs_T["Z"]
    E_avgs[T_ind, h_ind] = avgs_T["E"]
    Esq_avgs[T_ind, h_ind] = avgs_T["Esq"]
    M_avgs[T_ind, h_ind] = avgs_T["M"]
    Msq_avgs[T_ind, h_ind] = avgs_T["Msq"]
    N_tot_avgs[T_ind, h_ind] = avgs_T["N"]
    

    return Dict("Z" => Zs, "E" => E_avgs, "Esq" => Esq_avgs, "M" => M_avgs, "Msq" => Msq_avgs, "N" => N_tot_avgs)
end
"""

#=
Main.include("Hubbard.jl")
bonds=[[1,2]]
N=2; NTOP=0
sectors_info = sectors_info_gen(N=N)
quantities = E_Quants(N=N, U=10.0, t=1.0, sectors_info=sectors_info, bonds=bonds)
#printing(quantities; Quant_names = quant_names, name = "test", NTOP=NTOP)
thermal_avg(;T=100,μ=5.0, name = "test",N=N, NTOP=NTOP)=#
