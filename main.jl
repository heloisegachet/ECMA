include("dualisation.jl")
include("PL_statique.jl")

function main()
    #PL_statique("data/10_ulysses_3.tsp")
    dualisation("data/10_ulysses_3.tsp")
end

main()