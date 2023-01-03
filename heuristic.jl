
function sort_nodes(W)
	# trie les sommets par poids décroissant
	sorted_nodes = sortperm(W, rev=true)
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

function heuristic(G, W, K, B)
	# retourne une solution realisable ou affiche une erreur
	sorted_nodes = sort_nodes(W)
	partition = Vector{Vector{Int64}}()
	for k in 1:K
		push!(partition, Vector{Int64}())
	end
	for i in sorted_nodes
		for k in 1:K
			if allowed(i,partition[k], W, B)
				push!(partition[k],i)
				break
			else
				if k == K
					println("impossible de placer ", i)
					return
				end
			end
		end
		
	end
	return partition
end

# test

G = [0 1 2 ; 1 0 3 ; 2 3 0]
W = [3,1,2]
K = 2
B = 3
println(heuristic(G,W,K,B))