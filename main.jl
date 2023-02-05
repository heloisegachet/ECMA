include("dualisation.jl")
include("PL_statique.jl")
include("plan_coupants.jl")
include("branch-and-cut.jl")
include("heuristic.jl")

function main()
    time_lim = 30
    folder = "data large/"
    PL_statique("data small/26_eil_6.tsp", 10*60)
    #foreach(readdir(folder)) do file
    #    println(file)
    #    println("statique")
    #    PL_statique(string(folder,file), time_lim)
    #    println("dualisation")
    #    dualisation(string(folder,file), time_lim)
    #    println("plan_coupants")
    #    plan_coupants(string(folder,file), time_lim)
    #    println("branch-and-cut")
    #    branch_and_cut(string(folder,file), time_lim)
    #    println("heuristic")
    #    heuristic(string(folder,file), time_lim)
    #end
    #foreach(readdir(folder)) do file
    #    println(file)
    #    println("statique")
    #    relax_statique = PL_statique(string(folder,file), time_lim)
    #    println("dualisation")
    #    dualisation(string(folder,file), time_lim, relax=relax_statique)
    #    println("plan_coupants")
    #    plan_coupants(string(folder,file), time_lim, relax=relax_statique)
    #    println("branch-and-cut")
    #    branch_and_cut(string(folder,file), time_lim, relax=relax_statique)
    #    println("heuristic")
    #    heuristic(string(folder,file), time_lim)
    #end
    #foreach(readdir(folder)) do file
    #    println(file)
    #    println("statique")
    #    PL_statique(string(folder,file), time_lim)
    #    println("heuristic")
    #    sol, value  = heuristic(string(folder,file), time_lim)
    #    println("plan_coupants")
    #    plan_coupants(string(folder,file), time_lim, sol_initiale=sol, val_initiale=value)
    #    println("branch-and-cut")
    #    branch_and_cut(string(folder,file), time_lim, sol_initiale=sol, val_initiale=value)
    #end

    #foreach(readdir(folder)) do file
    #    println(file)
    #    println("statique")
    #    PL_statique(string(folder,file), time_lim, cut_one=true)
    #    println("dualisation")
    #    dualisation(string(folder,file), time_lim, cut_one=true)
    #    println("plan_coupants")
    #    plan_coupants(string(folder,file), time_lim, cut_one=true)
    #    println("branch-and-cut")
    #    branch_and_cut(string(folder,file), time_lim, cut_one=true)
    #end

    #file = "10_ulysses_6.tsp"
    #PL_statique(string(folder,file), time_lim)
    #println("dualisation")
    #dualisation(string(folder,file), time_lim)
    #println("plan_coupants")
    #plan_coupants(string(folder,file), time_lim)
    #println("branch-and-cut")
    #branch_and_cut(string(folder,file), time_lim)
end
main()

