function added_value(i, partie, l)
	# augmentation du coût de la solution en ajoutant i à la partie
	v = 0
	for j in partie
		v = v + l[i,j]
	end
	return v
end


function heuristic(filepath)
	include(filepath)

	# liste des pires poids des sommets
	w_nodes = Vector{Float64}(zeros(n))
	for i in 1:n
		w_nodes[i] = w_v[i]*(1 + W_v[i])
	end

	# matrice des pires longueurs des arêtes
	l = zeros(Float64, n, n)
	for i in 1:n
		for j in 1:n
			if j!=i
				l[i,j] = sqrt((coordinates[i,1] - coordinates[j,1])^2 +(coordinates[i,2] - coordinates[j,2])^2) +3*(lh[i]+lh[j])
			end
		end
	end

	# constante de majoration (somme des longueurs)
	M = 0
	for i in 1:n
		for j in i+1:n
			M = M + l[i,j]
		end
	end

	# ordre des sommets par pire poids décroissant
	sorted_nodes = sortperm(w_nodes, rev=true)

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

	start = time()

	# heuristique goutonne
	for i in sorted_nodes
		# on place i au meilleur endroit possible
		best_k = 0
		min_added_value = M
		for k in 1:K
			if w_group[k] + w_nodes[i] > B
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
			partition = []
			v_sol = 0
			break
		else
			# placer i dans la partie best_k et maj des valeurs
			push!(partition[best_k], i)
			assignment[i] = best_k
			w_group[best_k] = w_group[best_k] + w_nodes[i]
			v_group[best_k] = v_group[best_k] + min_added_value
			v_sol = v_sol + min_added_value
		end
	end

	stop = time()

	# ecrire la solution dans un fichier
	fout = open("./solution_heuristique.txt", "a")
	if v_sol == 0
		println(fout, "file "*filepath*" with n = ", n, "erreur", "solving time = ",stop - start, "s")
	else
		println(fout, "file "*filepath*" with n = ", n, " solution of obj value ", v_sol, " solving time = ",stop - start, "s")
		println(fout, "solution : ")
		println(fout, partition)
	end
	close(fout)
end