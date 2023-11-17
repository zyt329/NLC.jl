"""
    prints the Quantities to a .dat file for save in a given folder for each cluster. Saved name is name_Ntop
    Input:
        Quantities::Array{Array{Float64,1},1} : An array of Arrays of quantities
            example: sample = [Float64[], Float64[], Float64[], Float64[], Float64[]]
                for i in 1:5
                    push!(sample[i], i)
                    push!(sample[i], 2i)
                end
        header::Array{Any,2} : An array of names of quantities
            example: ['E' "Esq" 'M' "Msq" "N_tot"]
        name::string : name of the saved files
        NTOP::Int64 : topological number of the cluster whose quantities we're saving.
"""
function printing(Quantities; name::String, NTOP::Int64, N::Int64, folder_path::String="./")

    header = ['E' "Esq" 'M' "Msq" "N_tot"]

    open(
        joinpath(folder_path, name * "_order$(N)_topo$(NTOP)" * ".dat"),
        "w",
    ) do io
        writedlm(io, header)
        for (i, Quant) in enumerate(Quantities[1])
            writedlm(io, [Quantities[n][i] for n in 1:length(Quantities)]')
        end
    end

end

# printing(Quantities; name = "test_print2txt", NTOP=1)

"""
    Reading in quantities of a cluster of order N and topological number NTOP, from the files with the name passed in.
"""
function reading_quantities(; name::String, NTOP::Int64, N::Int64, folder_path::String="./")

    global quantities

    # open folder and read
    open(
        joinpath(folder_path, name * "_order$(N)_topo$(NTOP)" * ".dat"),
        "r",
    ) do io
        quantities = readdlm(io, skipstart=1) # skip the first line
    end

    return quantities
end

"""
    Make folder that have non-repeating index by the end of its name. Give the folder an "sID" index (simulation Index) automatically if not specified.

    Returns the full path to the folder being created.
"""
function make_indexed_folder(; folder_prefix, folder_path=".", make_new=true, sID::Int=0)
    # initialize data folder names
    folder_name = @sprintf "%s-%d" folder_prefix sID
    folder_full_path = joinpath(folder_path, folder_name)

    if make_new
        # if null data folder id given, determine data name and id
        if sID == 0
            while isdir(folder_full_path) || sID == 0
                sID += 1
                folder_name = @sprintf "%s-%d" folder_prefix sID
                folder_full_path = joinpath(folder_path, folder_name)
            end
        end

        # make the folder
        mkpath(folder_full_path)

        return folder_full_path
    else
        # if null data folder id given, determine data name and id
        if sID == 0
            while isdir(folder_full_path) || sID == 0
                folder_name = @sprintf "%s-%d" folder_prefix sID
                folder_full_path = joinpath(folder_path, folder_name)
                sID += 1
            end
        end

        return folder_full_path
    end
end

"""
    Check how many sites there are in a cluster given bonds.
"""
function num_sites(bonds::Vector{Vector{Int}})
    sites = Set()
    for bond in bonds
        for site in bond
            push!(sites, site)
        end
    end
    return length(sites)
end


