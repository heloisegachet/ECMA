using JuMP
using CPLEX


function added_value(i, partie, l)
	# augmentation du coût de la solution en ajoutant i à la partie
	v = 0
	for j in partie
		v = v + l[i,j]
	end
	return v
end

function value_sol(partition, l)
	
	# Create the model
	m = Model(CPLEX.Optimizer)
	### Variables
	@variable(m, delta[i in 1:n, j in 1:n]>=0)
	### Constraints
	@constraint(m, [i in 1:n, j in 1:n;i!=j], delta[i,j]<=3)
	@constraint(m, sum(delta[i,j] for i in 1:n, j in 1:n if i!=j) <=L)
	### Objective
	obj = 0
	for index_k in 1:length(partition)
		if (length(partition[index_k]) >= 2)
			obj += sum(l[i,j] - 3*(lh[i]+lh[j]) + delta[i,j]*(lh[i]+lh[j]) for i in partition[index_k], j in partition[index_k] if i!=j) 
		end
	end
	
	@objective(m, Max, obj)
	### Solve the problem
	start = time()
	optimize!(m)
	stop = time()
	println(" solution of obj value ", JuMP.objective_value(m),"solving time = ",stop - start, "s")
	#for i in 1:n
	#	for j in 1:n
	#		println(i," ",j," l ", l[i,j]-3*(lh[i]+lh[j]), " delta ", value(delta[i,j]))
	#	end
	#end
	return JuMP.objective_value(m)
end


function heuristic(filepath)
	include(filepath)

	# liste des pires poids des sommets
	w_nodes = Vector{Float64}(zeros(n))
	for i in 1:n
		w_nodes[i] = w_v[i]*(1 + W_v[i])
	end
	# ordre des sommets par pire poids décroissant
	sorted_nodes = sortperm(w_nodes, rev=true)

	# matrice des pires longueurs des arêtes
	l = zeros(Float64, n, n)
	for i in 1:n
		for j in 1:n
			if j!=i
				l[i,j] = sqrt((coordinates[i,1] - coordinates[j,1])^2 +(coordinates[i,2] - coordinates[j,2])^2) +3*(lh[i]+lh[j])
			end
		end
	end

	# initialise la partition à vide
	partition = Vector{Vector{Int64}}()
	for k in 1:K
		push!(partition, Vector{Int64}())
	end
	# initialise le tableau indiquant pour chaque sommet sa partition (0 si non assigné)
	assignment = Vector{Int64}(zeros(n))
	# initialise le poids courant de chaque partie
	w_group = Vector{Float64}(zeros(K))
	# initialise la valeur courante de chaque partie
	v_group = Vector{Float64}(zeros(K))
	# initialise la valeur de la solution
	v_sol = 0
	# constante de majoration (somme des longueurs)
	M = 0
	for i in 1:n
		for j in i+1:n
			M = M + l[i,j]
		end
	end

	# heuristique goutonne
	for i in sorted_nodes
		# on place i au meilleur endroit possible
		best_k = 0
		min_added_value = M
		for k in 1:K
			if w_group[k] + w_nodes[i] > B
				continue
			else
				#println(partition[k],l)
				v = added_value(i, partition[k], l)
				if v < min_added_value
					min_added_value = v
					best_k = k
				end
			end
		end
		if best_k == 0
			println("erreur")
			partition = []
			v_sol = 0
			return partition, v_sol
		else
			# placer i dans la partie best_k et maj des valeurs
			push!(partition[best_k], i)
			assignment[i] = best_k
			w_group[best_k] = w_group[best_k] + w_nodes[i]
			v_group[best_k] += min_added_value
			v_sol += min_added_value
		end
	end
	v_sol = value_sol(partition, l)
	return partition, v_sol
end

filename = "data_small/10_ulysses_3.tsp"
P, v = heuristic(filename)
if v == 0
	println("erreur")
else
	println("partition = ", P)
	println("valeur = ", v)
end

# application à toutes les instances

#fout = open("output_heuristic.txt", "w")
#
#for file in readdir("./data")
#	println(fout, "instance " * file)
#	P, v = heuristic("./data/" * file)
#	if v == 0
#		println(fout, "erreur")
#	else
#		println(fout, "partition = ", P)
#		println(fout, "valeur = ", v)
#	end
#	println(fout, "")
#end
#
#close(fout)
#