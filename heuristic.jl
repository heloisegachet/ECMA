function added_value(i, partie, l)
	# calcule l'augmentation du coût de la solution en ajoutant i à la partie K
	v = 0
	for j in partie
		v = v + l[i,j]
	end
	return v
end


function solver_heuristic(filepath)
	include(filepath)

	# construit la liste des pires poids des sommets
	w_nodes = Vector{Float64}(zeros(n))
	for i in 1:n
		w_nodes[i] = w_v[i]+W_v[i]
	end

	# construit la matrice des pires longueurs des arêtes
	l = zeros(Float64, n, n)
	for i in 1:n
		for j in 1:n
			l[i,j] = sqrt((coordinates[i,1] - coordinates[j,1])^2 +(coordinates[i,2] - coordinates[j,2])^2) +3*(lh[i]+lh[j])
		end
	end
	
	# construit la liste des arêtes (i,j, lij) par pire longueur décroissante
	E = Vector{Tuple{Int64, Int64, Float64}}()
	for i in 1:n
		for j in i+1:n
			push!(E, (i,j,l[i,j]))
		end
	end
	sorted_E = sort(E, by = x -> x[3], rev=true)

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

	# heuristique gloutonne
	for e in sorted_E
		nodes = sort([(e[1], w_nodes[e[1]]), (e[2], w_nodes[e[2]])], by = x -> x[2], rev=true)
		for (i,p) in nodes
			if assignment[i]!=0
				continue
			else
				# on place i au meilleur endroit possible
				best_k = 0
				min_added_value = 100000000
				for k in 1:K
					if w_group[k]+w_nodes[i] > B
						continue
					else
						v = added_value(i, partition[k], l)
						if v < min_added_value
							min_added_value = v
							best_k = k
						end
					end
				end
				if best_k == 0
					println("erreur")
					return
				else
					# placer i dans la partie best_k et maj des valeurs
					push!(partition[best_k], i)
					assignment[i] = best_k
					w_group[best_k] = w_group[best_k] + w_nodes[i]
					v_group[best_k] = v_group[best_k] + min_added_value
					v_sol = v_sol + min_added_value
				end
			end
		end
	end
	println(partition)
	println(v_sol)
end

# calcul pour une instance

solver_heuristic("./data/10_ulysses_3.tsp")