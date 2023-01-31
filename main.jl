include("dualisation.jl")
include("PL_statique.jl")
include("plan_coupants.jl")
include("branch-and-cut.jl")
include("heuristic.jl")

function main()
    #PL_statique("data/10_ulysses_3.tsp")
    #dualisation("data/10_ulysses_3.tsp")
    #print("z star final ", plan_coupants("data/10_ulysses_3.tsp"))
    #branch_and_cut("data/10_ulysses_3.tsp")

    #foreach(readdir("data_small/")) do file
    #    PL_statique(string("data_small/",file))
    #end
    foreach(readdir("data_small/")) do file
        dualisation(string("data_small/",file))
    end
    #partition, value = heuristic("data/10_ulysses_3.tsp")
    #println(partition, value)

    #println("Branch and cut warm start")
    #branch_and_cut("data/22_ulysses_6.tsp", partition, value)
    #println("Branch and cut")
    #branch_and_cut("data/10_ulysses_6.tsp", [], 0)

    #dualisation("data/22_ulysses_6.tsp", partition)
    #dualisation("data/22_ulysses_6.tsp", [])
end

main()

