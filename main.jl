include("dualisation.jl")
include("PL_statique.jl")
include("plan_coupants.jl")

function main()
    #PL_statique("data/10_ulysses_3.tsp")
    dualisation("data/10_ulysses_3.tsp")
    #print("z star final ", plan_coupants("data/10_ulysses_3.tsp"))

    #foreach(readdir("data_small/")) do file
    #    PL_statique(string("data_small/",file))
    #end
    #foreach(readdir("data_small/")) do file
    #    dualisation(string("data_small/",file))
    #end
end

main()

