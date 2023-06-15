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

    # loop over all symmetry sectors of Sᶻ=m           
    # only loop over half of the sectors
    # other sectors related by flipping all spins (ΠᵢSˣᵢ)
    for m = 0:Int(floor(N))

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
function thermal_avg(; T::Real, J::Real, quantities, h::Real=0.0, g::Real)
    # passing in quantities
    E::Array{Float64,1} = quantities[1]
    Esq::Array{Float64,1} = quantities[2]
    M::Array{Float64,1} = quantities[3]
    Msq::Array{Float64,1} = quantities[4]
    N_tot::Array{Float64,1} = quantities[5]
    # calculate thermal average
    β = 1 / T
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
    return [Z E_avg Esq_avg M_avg Msq_avg N_tot_avg]
end

"""
    Read in quantities(eigen energies and corresponding average quantities of eigenstates), a range of temperatures at which we want to calculate the thermal average. Out put thermal average of quantities at the read-in temperatures.
"""
function thermal_avg_hT_loop(; Temps, J::Real, quantities, hs=0.0, g::Real=2.1)
    if hs == 0.0 # if "hs" is not passed in
        Zs = Float64[]
        E_avgs = Float64[]
        Esq_avgs = Float64[]
        M_avgs = Float64[]
        Msq_avgs = Float64[]
        N_tot_avgs = Float64[]
        # don't loop over h if hs is not passed in
        for T in Temps
            avgs_T = thermal_avg(; T=T, J=J, quantities=quantities, h=0.0, g=g)
            push!(Zs, avgs_T[1])
            push!(E_avgs, avgs_T[2])
            push!(Esq_avgs, avgs_T[3])
            push!(M_avgs, avgs_T[4])
            push!(Msq_avgs, avgs_T[5])
            push!(N_tot_avgs, avgs_T[6])
        end
    else # if "hs" is passed in
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
                avgs_T = thermal_avg(T=T, J=J, quantities=quantities, h=h, g=g)
                Zs[T_ind, h_ind] = avgs_T[1]
                E_avgs[T_ind, h_ind] = avgs_T[2]
                Esq_avgs[T_ind, h_ind] = avgs_T[3]
                M_avgs[T_ind, h_ind] = avgs_T[4]
                Msq_avgs[T_ind, h_ind] = avgs_T[5]
                N_tot_avgs[T_ind, h_ind] = avgs_T[6]
            end
        end
    end
    return [Zs, E_avgs, Esq_avgs, M_avgs, Msq_avgs, N_tot_avgs]
end

quant_names = ['E' "Esq" 'M' "Msq" "N_tot"]
#=
Main.include("Hubbard.jl")
bonds=[[1,2]]
N=2; NTOP=0
sectors_info = sectors_info_gen(N=N)
quantities = E_Quants(N=N, U=10.0, t=1.0, sectors_info=sectors_info, bonds=bonds)
#printing(quantities; Quant_names = quant_names, name = "test", NTOP=NTOP)
thermal_avg(;T=100,μ=5.0, name = "test",N=N, NTOP=NTOP)=#
