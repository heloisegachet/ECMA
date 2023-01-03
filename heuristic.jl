
function sort_nodes(poids)
	# trie les sommets par poids décroissant
	sorted_nodes = sortperm(poids, rev=true)
	return sorted_nodes
end

function allowed(i, partie, W, B)
	# verifie qu'on peut ajouter le sommet i dans la partie
	w_partie = 0
	for j in partie
		w_partie += W[j]
	end
	w = W[i]
	if w_partie + w <= B
		return true
	else
		return false
	end
end

function solver_heuristic(filepath)
	include(filepath)

	# construit la liste des poids dans le pire des cas
	poids = Vector{Float64}()
	for i in 1:n
		push!(poids, w_v[i]+W_v[i])
	end

	# appel l'heuristique sur ces donnees
	return heuristic(poids, K, B)
end

function heuristic(poids, K, B)
	# retourne une solution realisable ou affiche une erreur
	
	sorted_nodes = sort_nodes(poids)
	partition = Vector{Vector{Int64}}()
	for k in 1:K
		push!(partition, Vector{Int64}())
	end
	for i in sorted_nodes
		for k in 1:K
			if allowed(i,partition[k], poids, B)
				push!(partition[k],i)
				break
			else
				if k == K
					println("erreur")
					return
				end
			end
		end
		
	end
	return partition
end

# test

for file in readdir("./data")
	# Afficher le nom du fichier
	println(solver_heuristic("./data/"*file))
end
