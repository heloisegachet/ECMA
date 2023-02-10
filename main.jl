include("dualisation.jl")
include("PL_statique.jl")
include("plan_coupants.jl")
include("branch-and-cut.jl")
include("heuristic.jl")

function main()
    time_lim = 60
    folder = "data/"
    foreach(readdir(folder)) do file
        println(file)
        println("statique")
        #PL_statique(string(folder,file), time_lim)
        println("dualisation")
        res = dualisation(string(folder,file), time_lim, gap=1e-6)
        if !isnothing(res)
            sol, val, gap =res
        else
            gap = nothing
        end
        println("plan_coupants")
        #plan_coupants(string(folder,file), time_lim)
        println("branch-and-cut")
        #branch_and_cut(string(folder,file), time_lim)
        println("heuristic")
        if !isnothing(gap) && gap < 1e-8
            heuristic(string(folder,file), time_lim, val_min=val)
        else
            heuristic(string(folder,file), time_lim, val_min=nothing)
        end
    end

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
