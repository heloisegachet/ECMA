include("dualisation.jl")
include("PL_statique.jl")
include("plan_coupants.jl")
include("branch-and-cut.jl")
include("heuristic.jl")

function main()
    time_lim = 120
    folder = "data small/"
    foreach(readdir(folder)) do file
        println(file)
        println("statique")
        PL_statique(string(folder,file), time_lim)
        println("dualisation")
        dualisation(string(folder,file), time_lim)
        println("plan_coupants")
        plan_coupants(string(folder,file), time_lim)
        println("branch-and-cut")
        branch_and_cut(string(folder,file), time_lim)
        println("heuristic")
        heuristic(string(folder,file), time_lim)
    end
end
main()

