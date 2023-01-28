include("dualisation.jl")
include("PL_statique.jl")
include("plan_coupants.jl")
include("heuristic.jl")
include("k-center.jl")
include("k-meanspp.jl")

function main()
	#PL_statique("data/10_ulysses_3.tsp")
	#dualisation("data/10_ulysses_3.tsp")
	#print("z star final ", plan_coupants("data/10_ulysses_3.tsp"))
	#heuristic("./data/10_ulysses_3.tsp")
end

main()

