using Clustering
using Random


function value(partie, C)
	v_p = 0
	m = length(partie)
	for i in 1:m
		for j in i+1:m
			v_p = v_p + C[i,j]
		end
	end
	return v_p
end


function added_value(i, partie, C)
	# augmentation du coût de la solution en ajoutant i à la partie
	v = 0
	for j in partie
		v = v + C[i,j]
	end
	return v
end


function isFeasible(partition, w_nodes, B)
	for p in partition
		w_p = 0
		for i in p
			w_p = w_p + w_nodes[i]
		end
		if w_p > B
			return false
		end
	end
	return true
end


function k_means_pp(filepath)
	include(filepath)

	# liste des pires poids des sommets
	w_nodes = Vector{Float64}(zeros(n))
	for i in 1:n
		w_nodes[i] = w_v[i]*(1 + W_v[i])
	end

	# formatage des données
	X = transpose(coordinates)
	# matrice des coûs = pires longueurs
	C = zeros(Float64, n, n)
	for i in 1:n
		for j in 1:n
			if j!=i
				C[i,j] = sqrt((coordinates[i,1] - coordinates[j,1])^2 +(coordinates[i,2] - coordinates[j,2])^2) +3*(lh[i]+lh[j])
			end
		end
	end

	# partitionnement
	R = kmeans(X, K; maxiter=200, display=:none, init=kmpp_by_costs(C, K))
	@assert nclusters(R) == K
	a = assignments(R)

	# contruit la partition associée
	partition = Vector{Vector{Int64}}()
	for k in 1:K
		push!(partition, Vector{Int64}())
	end
	for i in 1:n
		push!(partition[a[i]], i)
	end
	
	# réparation
	feasible = isFeasible(partition, w_nodes, B)

	if feasible == false

		sorted_nodes = sortperm(w_nodes, rev=true)
		w_group = Vector{Float64}(zeros(K))
		p_fixed = Vector{Vector{Int64}}()
		for k in 1:K
			push!(p_fixed, Vector{Int64}())
		end
		v_group = Vector{Float64}(zeros(K))
		v_sol = 0
		M = 0
		for i in 1:n
			for j in i+1:n
				M = M + C[i,j]
			end
		end

		for i in sorted_nodes
			if w_group[a[i]]+w_nodes[i]<B
				push!(p_fixed[a[i]], i)
				w_group[a[i]] = w_group[a[i]]+w_nodes[i]
				add_v = added_value(i, p_fixed[a[i]], C)
				v_group[a[i]] = v_group[a[i]] + add_v
				v_sol = v_sol + add_v
			else
				# on met i là où on peut, tout de suite
				best_k = 0
				min_added_value = M
				for k in 1:K
					if w_group[k] + w_nodes[i] > B
						continue
					else
						v = added_value(i, p_fixed[k], C)
						if v < min_added_value
							min_added_value = v
							best_k = k
						end
					end
				end
				if best_k == 0
					println("erreur")
					p_fixed = []
					v_sol = 0
					break
				else
					# placer i dans la partie best_k et maj des valeurs
					push!(p_fixed[best_k], i)
					w_group[best_k] = w_group[best_k] + w_nodes[i]
					v_group[best_k] = v_group[best_k] + min_added_value
					v_sol = v_sol + min_added_value
				end
			end
		end

		return p_fixed, v_sol
	end
	v_group = Vector{Float64}(zeros(K))
	v_sol = 0
	for k in 1:K
		v_group[k] = value(partition[k], C)
		v_sol = v_sol + v_group[k]
	end
	return partition, v_sol
end


function test()
	fout = open("output_k_means_pp.txt", "w")
	for file in readdir("./data")
		println(fout, "instance " * file)
		seed = 1
		Random.seed!(seed)
		p, v = k_means_pp("./data/"*file)
		while(v == 0 && seed < 3)
			seed = seed +1
			Random.seed!(seed)
			p, v = k_means_pp("./data/"*file)
		end
		println(fout, "partition = ", p)
		println(fout, "valeur = ", v)
		
	end
	close(fout)
end