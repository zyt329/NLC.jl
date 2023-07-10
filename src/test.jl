using JLD2

jldopen("/nfs/home/zyt329/Research/xxz/NLC.jl/src/testf" * ".jld2",
    "a+", compress=true) do file
    for N in 1:10^2

        file["E"*string(N)] = rand(100, 100)
        file["M"*string(N)] = rand(100, 100)

    end
end

file = jldopen("/nfs/home/zyt329/Research/xxz/NLC.jl/src/testf" * ".jld2", "r")

for N in 1:100
    println("N=$N")
    E = file["E"*string(N)]
end

