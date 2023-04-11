using LinearAlgebra
using SparseArrays
using Arpack

"""
    Check for the state, how many spins points up.
    Input: state: the state to check.
           N: the number of sites in the cluster.

return the number of up spins in the state.
"""
function chk_m(state::Int64,N::Int64)
    state_binary = digits!(zeros(Int64, 64), state, base = 2)[1:N]
    return sum(state_binary)
end

"""
    Loop over all states, check which sector the state belongs.
return a dictionary:
    :states::Array{Array{Int,1},1} : the m+1 th entry gives an array of states in the m sector
    :state_tot::Array{Int,1} : the m+1 th entry gives the total number of states in the m sector.
    :state_num::Dict{Int64, Int64} : give the state as key, returns the numbering of the state in its m sector.
"""
function sectors_info_gen(;N::Int64)
    states::Array{Array{Int,1},1} = Array{Int,1}[[] for i in 0:N]
    state_tot::Array{Int,1} = Int[0 for i in 0:N]
    state_num::Dict{Int64, Int64} = Dict{Int64, Int64}()
    for state in 0:(2^N-1)
        m = chk_m(state, N)
        push!(states[m+1], state)
        state_tot[m+1] += 1
        state_num[state] = state_tot[m+1]
    end
    @assert(sum(state_tot) == 2^N, "total number of states is not 2^N")
    return Dict{Symbol, Any}(:states => states, :state_tot => state_tot, :state_num => state_num)
end

using SparseArrays

function update_val(row_inds, col_inds, vals;row_ind, col_ind, val)
    push!(row_inds, row_ind)
    push!(col_inds, col_ind)
    push!(vals, val)
end

function H_sector(;J_xy::Float64, J_z::Float64, N::Int64, m::Int64, sectors_info::Dict{Symbol,Any}, bonds)
    row_inds = Int64[]
    col_inds = Int64[]
    vals = Float64[]
    states = sectors_info[:states][m+1]
    state_tot = sectors_info[:state_tot][m+1]
    state_num = sectors_info[:state_num]
    for state in states #loop over all states in the sector
        state_binary = digits!(zeros(Int64, 64), state, base = 2)
        for bond in bonds #bond=[s1, s2], where s1,s2 are the two sites of the bond
            s1 = bond[1]; s2 = bond[2]
            if state_binary[s1] == state_binary[s2]
                update_val(row_inds, col_inds, vals, row_ind = state_num[state], col_ind = state_num[state], val = (1/4) * J_z)
            else
                update_val(row_inds, col_inds, vals, row_ind = state_num[state], col_ind = state_num[state], val = -(1/4) * J_z)
                flipped_state = state ⊻ (1<<(s1-1))
                flipped_state = flipped_state ⊻ (1<<(s2-1))
                update_val(row_inds, col_inds, vals, row_ind = state_num[state], col_ind = state_num[flipped_state], val = (1/2) * J_xy)
            end
        end
    end
    return sparse(row_inds, col_inds, vals, state_tot, state_tot, +)
end
