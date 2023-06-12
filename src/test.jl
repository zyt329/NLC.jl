using JLD

save(
    "test.jld",
    "dict",
    Dict("1" => 1)
)