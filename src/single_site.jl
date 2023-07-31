# calculate properties for single-site xxz model
# Should import quantities form python script
include("./consts.jl")

function single_site_quantities_Heisbg(;Ts, hs=0.0, g::Real)
    if hs==0.0 # if hs is not passed in
        Z = Float64[]
        E = Float64[]
        Esq = Float64[]
        M = Float64[]
        Msq = Float64[]
        N_tot = Float64[]
        for T in Ts
            β = 1/T
            Zt = exp(0)+exp(0)
            push!(Z, Zt)
            push!(E,( 0) / Zt)
            push!(Esq,(0) / Zt)
            push!(M,0)
            push!(Msq , 2*(1/4*g^2*mu_B^2)/Zt)
            push!(N_tot , 2*(1)/Zt)
        end
    else # if hs is passed in
        NT = length(Ts)
        Nh = length(hs)
        Z = Array{Float64}(undef, NT, Nh)
        E = Array{Float64}(undef, NT, Nh)
        Esq = Array{Float64}(undef, NT, Nh)
        M = Array{Float64}(undef, NT, Nh)
        Msq = Array{Float64}(undef, NT, Nh)
        N_tot = Array{Float64}(undef, NT, Nh)
        for (h_ind, h) in enumerate(hs)
            for (T_ind, T) in enumerate(Ts)
                β = 1/T
                Zt = exp(β*(h*g*mu_B/2/k_B))+exp(-β*(h*g*mu_B/2/k_B))
                Z[T_ind, h_ind] = Zt
                E[T_ind, h_ind] = ((-h*g*mu_B/2)*exp(β*(h*g*mu_B/2/k_B))+(h*g*mu_B/2)*exp(-β*(h*g*mu_B/2/k_B)))/ Zt
                Esq[T_ind, h_ind] = ((-h*g*mu_B/2)^2*exp(β*(h*g*mu_B/2/k_B))+(h*g*mu_B/2)^2*exp(-β*(h*g*mu_B/2/k_B)))/ Zt
                M[T_ind, h_ind] = ((-g*mu_B/2)*exp(β*(h*g*mu_B/2/k_B))+(g*mu_B/2)*exp(-β*(h*g*mu_B/2/k_B)))/ Zt
                Msq[T_ind, h_ind] = ((-g*mu_B/2)^2*exp(β*(h*g*mu_B/2/k_B))+(g*mu_B/2)^2*exp(-β*(h*g*mu_B/2/k_B)))/ Zt
                N_tot[T_ind, h_ind] = 1
            end
        end
    end
    return [Z, E,Esq, M, Msq, N_tot]
end

"""
Should be the same as Heisbg model.
"""
function single_site_quantities_xxz(;Ts, hs=0.0, g::Real)
    if hs==0.0 # if hs is not passed in
        Z = Float64[]
        E = Float64[]
        Esq = Float64[]
        M = Float64[]
        Msq = Float64[]
        N_tot = Float64[]
        for T in Ts
            β = 1/T
            Zt = exp(0)+exp(0)
            push!(Z, Zt)
            push!(E,( 0) / Zt)
            push!(Esq,(0) / Zt)
            push!(M,0)
            push!(Msq , 2*(1/4*g^2*mu_B^2)/Zt)
            push!(N_tot , 2*(1)/Zt)
        end
    else # if hs is passed in
        NT = length(Ts)
        Nh = length(hs)
        Z = Array{Float64}(undef, NT, Nh)
        E = Array{Float64}(undef, NT, Nh)
        Esq = Array{Float64}(undef, NT, Nh)
        M = Array{Float64}(undef, NT, Nh)
        Msq = Array{Float64}(undef, NT, Nh)
        N_tot = Array{Float64}(undef, NT, Nh)
        for (h_ind, h) in enumerate(hs)
            for (T_ind, T) in enumerate(Ts)
                β = 1/T
                Zt = exp(β*(h*g*mu_B/2/k_B))+exp(-β*(h*g*mu_B/2/k_B))
                Z[T_ind, h_ind] = Zt
                E[T_ind, h_ind] = ((-h*g*mu_B/2)*exp(β*(h*g*mu_B/2/k_B))+(h*g*mu_B/2)*exp(-β*(h*g*mu_B/2/k_B)))/ Zt
                Esq[T_ind, h_ind] = ((-h*g*mu_B/2)^2*exp(β*(h*g*mu_B/2/k_B))+(h*g*mu_B/2)^2*exp(-β*(h*g*mu_B/2/k_B)))/ Zt
                M[T_ind, h_ind] = ((-g*mu_B/2)*exp(β*(h*g*mu_B/2/k_B))+(g*mu_B/2)*exp(-β*(h*g*mu_B/2/k_B)))/ Zt
                Msq[T_ind, h_ind] = ((-g*mu_B/2)^2*exp(β*(h*g*mu_B/2/k_B))+(g*mu_B/2)^2*exp(-β*(h*g*mu_B/2/k_B)))/ Zt
                N_tot[T_ind, h_ind] = 1
            end
        end
    end
    return Dict("Z" => Z, "E" => E, "Esq" => Esq, "M" => M, "Msq" => Msq, "N" => N_tot)
    #[Z, E,Esq, M, Msq, N_tot]
end
#result = single_site_quantities(Ts=Temps,hs=[-2.0],g=2.1)
#println(result[4]*N_A)
#=
using Plots

E = single_site_quantities(Ts=range(0.1,10,length=100),μ=0.0,U=10.0)[1]

plot(range(0.1,10,length=100), E)=#
