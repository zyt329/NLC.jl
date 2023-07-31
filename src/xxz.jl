using LinearAlgebra
using SparseArrays
using Arpack

"""
    Check for the state, how many spins points up.
    Input: state: the state to check.
           N: the number of sites in the cluster.

return the number of up spins in the state.
"""
function chk_m(state::Int64, N::Int64)
    state_binary = digits!(zeros(Int64, 64), state, base=2)[1:N]
    return sum(state_binary)
end

"""
    Loop over all states, check which sector the state belongs.
return a dictionary:
    :states::Array{Array{Int,1},1} : the m+1 th entry gives an array of states in the m sector
    :state_tot::Array{Int,1} : the m+1 th entry gives the total number of states in the m sector.
    :state_num::Dict{Int64, Int64} : give the state as key, returns the numbering of the state in its m sector.
"""
function sectors_info_gen(; N::Int64)
    states::Array{Array{Int,1},1} = Array{Int,1}[[] for i in 0:N]
    state_tot::Array{Int,1} = Int[0 for i in 0:N]
    state_num::Dict{Int64,Int64} = Dict{Int64,Int64}()
    for state in 0:(2^N-1)
        m = chk_m(state, N)
        push!(states[m+1], state)
        state_tot[m+1] += 1
        state_num[state] = state_tot[m+1]
    end
    @assert(sum(state_tot) == 2^N, "total number of states is not 2^N")
    return Dict{Symbol,Any}(:states => states, :state_tot => state_tot, :state_num => state_num)
end

"""
    Generate sector info for Sᶻ=0 (m=0) sector when N is even, using ΠᵢSˣᵢ symmetry.
"""
function m0_sectors_info_gen(; sectors_info::Dict{Symbol,Any}, N::Int64)

    # assert that N is even
    @assert(iseven(N), "N must be even to have a total z-spin 0 eigen state!")

    # load in m=N/2 sector info
    states = sectors_info[:states][Int(N / 2)+1]


    # holder of even/odd states
    even_states = Int64[]
    odd_states = Int64[]

    even_state_counter = 0
    odd_state_counter = 0

    even_state_num = Dict{Int64,Int64}()
    odd_state_num = Dict{Int64,Int64}()

    # loop over all m=N/2 sector's states
    for state in states

        # find the all-spin flipped state 
        # only flip the first N digits
        all_spin_flipped_state = 2^N - 1 - state

        # put states into the correct container
        # larger-valued state reps even state, smaller-valued state reps odd state
        if state > all_spin_flipped_state
            push!(even_states, state)
            #push!(odd_states, all_spin_flipped_state)

            even_state_counter += 1

            even_state_num[state] = even_state_counter
            #odd_state_num[all_spin_flipped_state] = state_tot
        else
            #push!(even_states, all_spin_flipped_state)
            push!(odd_states, state)

            odd_state_counter += 1

            #even_state_num[all_spin_flipped_state] = state_tot
            odd_state_num[state] = odd_state_counter
        end
    end
    @assert(even_state_counter + odd_state_counter == binomial(N, Int(N / 2)), "total number of states is not N choose N/2")

    return Dict{Symbol,Any}(:even_states => even_states, :odd_states => odd_states, :state_tot => odd_state_counter, :even_state_num => even_state_num, :odd_state_num => odd_state_num)
end


using SparseArrays

function update_val(row_inds, col_inds, vals; row_ind, col_ind, val)
    push!(row_inds, row_ind)
    push!(col_inds, col_ind)
    push!(vals, val)
end

function H_sector(; J_xy::Float64, J_z::Float64, N::Int64, m::Int64, sectors_info::Dict{Symbol,Any}, bonds)
    row_inds = Int64[]
    col_inds = Int64[]
    vals = Float64[]
    states = sectors_info[:states][m+1]
    state_tot = sectors_info[:state_tot][m+1]
    state_num = sectors_info[:state_num]
    for state in states #loop over all states in the sector
        state_binary = digits!(zeros(Int64, 64), state, base=2)
        for bond in bonds #bond=[s1, s2], where s1,s2 are the two sites of the bond
            s1 = bond[1]
            s2 = bond[2]
            if state_binary[s1] == state_binary[s2]
                update_val(row_inds, col_inds, vals, row_ind=state_num[state], col_ind=state_num[state], val=(1 / 4) * J_z)
            else
                update_val(row_inds, col_inds, vals, row_ind=state_num[state], col_ind=state_num[state], val=-(1 / 4) * J_z)
                flipped_state = state ⊻ (1 << (s1 - 1))
                flipped_state = flipped_state ⊻ (1 << (s2 - 1))
                update_val(row_inds, col_inds, vals, row_ind=state_num[state], col_ind=state_num[flipped_state], val=(1 / 2) * J_xy)
            end
        end
    end
    return sparse(row_inds, col_inds, vals, state_tot, state_tot, +)
end

"""
    Generate Hamiltonian for Sᶻ=0 (m=0) sector when N is even, using ΠᵢSˣᵢ symmetry.
"""
function H_m0_sector(; J_xy::Float64, J_z::Float64, N::Int64, m0_sectors_info::Dict{Symbol,Any}, bonds)

    # read in sector info
    even_states = m0_sectors_info[:even_states]
    odd_states = m0_sectors_info[:odd_states]
    even_state_num = m0_sectors_info[:even_state_num]
    odd_state_num = m0_sectors_info[:odd_state_num]

    state_tot = m0_sectors_info[:state_tot]

    # stores Hamiltonian info of even sectors
    row_inds_even = Int64[]
    col_inds_even = Int64[]
    vals_even = Float64[]

    # Hamiltonian for even sector
    for state in even_states #loop over all states in the sector

        state_binary = digits!(zeros(Int64, 64), state, base=2)
        for bond in bonds #bond=[s1, s2], where s1,s2 are the two sites of the bond
            s1 = bond[1]
            s2 = bond[2]
            if state_binary[s1] == state_binary[s2]
                # diagonal element
                update_val(row_inds_even, col_inds_even, vals_even, row_ind=even_state_num[state], col_ind=even_state_num[state], val=(1 / 4) * J_z)
            else
                # ============ diagonal element =============
                update_val(row_inds_even, col_inds_even, vals_even, row_ind=even_state_num[state], col_ind=even_state_num[state], val=-(1 / 4) * J_z)
                # ===========================================

                # ============ off diagonal element ==============
                flipped_state = state ⊻ (1 << (s1 - 1))
                flipped_state = flipped_state ⊻ (1 << (s2 - 1))

                # flip all the spins of the "flipped_state"
                reflected_flipped_state = 2^N - 1 - flipped_state

                if flipped_state > reflected_flipped_state # 'flipped_state' in the even sector
                    update_val(row_inds_even, col_inds_even, vals_even, row_ind=even_state_num[state], col_ind=even_state_num[flipped_state], val=(1 / 2) * J_xy)
                else # 'reflected_flipped_state' in the even sector
                    update_val(row_inds_even, col_inds_even, vals_even, row_ind=even_state_num[state], col_ind=even_state_num[reflected_flipped_state], val=(1 / 2) * J_xy)
                end
                # =================================================
            end
        end
    end

    # stores Hamiltonian info of odd sectors
    row_inds_odd = Int64[]
    col_inds_odd = Int64[]
    vals_odd = Float64[]

    # Hamiltonian for odd sector
    for state in odd_states #loop over all states in the sector

        state_binary = digits!(zeros(Int64, 64), state, base=2)
        for bond in bonds #bond=[s1, s2], where s1,s2 are the two sites of the bond
            s1 = bond[1]
            s2 = bond[2]
            if state_binary[s1] == state_binary[s2]
                # diagonal element
                update_val(row_inds_odd, col_inds_odd, vals_odd, row_ind=odd_state_num[state], col_ind=odd_state_num[state], val=(1 / 4) * J_z)
            else
                # ============ diagonal element =============
                update_val(row_inds_odd, col_inds_odd, vals_odd, row_ind=odd_state_num[state], col_ind=odd_state_num[state], val=-(1 / 4) * J_z)
                # ===========================================

                # ============ off diagonal element ==============
                flipped_state = state ⊻ (1 << (s1 - 1))
                flipped_state = flipped_state ⊻ (1 << (s2 - 1))

                # flip all the spins of the "flipped_state"
                reflected_flipped_state = 2^N - 1 - flipped_state

                if flipped_state < reflected_flipped_state # 'flipped_state' is basis state in the odd sector
                    update_val(row_inds_odd, col_inds_odd, vals_odd, row_ind=odd_state_num[state], col_ind=odd_state_num[flipped_state], val=(1 / 2) * J_xy)
                else # 'reflected_flipped_state' is basis state in the odd sector
                    update_val(row_inds_odd, col_inds_odd, vals_odd, row_ind=odd_state_num[state], col_ind=odd_state_num[reflected_flipped_state], val=-(1 / 2) * J_xy)
                end
                # =================================================
            end
        end
    end

    # return even_H, odd_H
    return sparse(row_inds_even, col_inds_even, vals_even, state_tot, state_tot, +), sparse(row_inds_odd, col_inds_odd, vals_odd, state_tot, state_tot, +)
end