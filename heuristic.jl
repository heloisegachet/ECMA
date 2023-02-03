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

function value_sol(partition)
	# Create the model
	m = Model(CPLEX.Optimizer)
	set_silent(m)
	### Variables
	@variable(m, delta[i in 1:n, j in 1:n]>=0)
	### Constraints
	@constraint(m, [i in 1:n, j in 1:n;i!=j], delta[i,j]<=3)
	@constraint(m, sum(delta[i,j] for i in 1:n, j in 1:n if i!=j) <=L)
	### Objective
	@objective(m, Max, sum(sum(l_statique[i,j] + delta[i,j]*(lh[i]+lh[j]) for i in partition[index_k], j in partition[index_k] if i<j) for index_k in 1:size(partition)[1] if (length(partition[index_k]) >= 2)))
	### Solve the problem
	optimize!(m)
	return objective_value(m)
end

function poids_partie(partie, time_lim)
	if time_lim < 1
		return
	end
	# Create the model
	m = Model(CPLEX.Optimizer)
	set_silent(m)
	set_time_limit_sec(m, time_lim)
	### Variables
	@variable(m, delta[v in partie]>=0)
	### Constraints
	@constraint(m, [v in partie], delta[v]<=W_v[v])
	@constraint(m, sum(delta[v] for v in partie) <=W)
	@objective(m, Max, sum(w_v[v]*(1+delta[v]) for v in partie))
	### Solve the problem
	optimize!(m)
	if(has_values(m))
		return objective_value(m)
	end
	return
end

function heuristic(filename, time_lim)

	include(filename)

	start = time()

	# liste des pires poids des sommets
	w_nodes = Vector{Float64}(zeros(n))
	for i in 1:n
		w_nodes[i] = w_v[i]*(1 + W_v[i])
	end

	# matrice des pires longueurs des arêtes
	global l_statique = zeros(Float64, n, n)
	l = zeros(Float64, n, n)

	for i in 1:n
		for j in 1:n
			if j!=i
				l_statique[i,j] = sqrt((coordinates[i,1] - coordinates[j,1])^2 +(coordinates[i,2] - coordinates[j,2])^2)
				l[i,j] = l_statique[i,j] +3*(lh[i]+lh[j])
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

	# heuristique goutonne
	for i in sorted_nodes
		# on place i au meilleur endroit possible
		best_k = 0
		min_added_value = M
		for k in 1:K
			new_partie = Vector{Int64}()
			for elem in partition[k]
				push!(new_partie,elem)
			end
			push!(new_partie, i)
			poids = poids_partie(new_partie, time_lim-time()+start)
			if !isnothing(poids) && poids > B 
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
			return
		else
			# placer i dans la partie best_k et maj des valeurs
			push!(partition[best_k], i)
			assignment[i] = best_k
			w_group[best_k] = w_group[best_k] + w_nodes[i]
			v_group[best_k] += min_added_value
		end
	end

	old_val = value_sol(partition)
	println("RECHERCHE LOCALE")
	best_partition, best_val = recherche_locale(partition, time_lim - time() + start)

	stop = time()

	write("heuristique", filename,stop-start, best_partition, string(old_val," / ",best_val), "None", "None")
	return best_partition, best_val
end


function recherche_locale(partition, time_lim)
	start = time()
	old_val = value_sol(partition)
	iter = 0
	while time() - start < time_lim
		k1 = rand(1:K)
		k2 = rand([1:k1-1;k1+1:K])
		v1 = rand([partition[k1];-1])
		v2 = rand(partition[k2])
		#println("sommets échanges ",v1," ",v2)
		new_partition = Vector{Vector{Int64}}()
		for k in 1:K
			push!(new_partition, Vector{Int64}())
			for i in partition[k]
				if (k!=k1 || i!=v1) && (k!=k2 || i!=v2)
					push!(new_partition[k], i)
				end
			end
			if k==k1
				push!(new_partition[k], v2)
			end
			if k==k2 && v1 != -1 
				push!(new_partition[k], v1)
			end
		end
		valid = true
		for partie in new_partition
			poids = poids_partie(partie, time_lim - time()+start)
			if isnothing(poids)
				return partition, old_val
			end
			if poids > B 
				valid = false
				break
			end
		end
		if valid
			#print("test valeur echange ")
			new_val = value_sol(new_partition)
			if old_val > new_val
				#println("better replace")
				partition = new_partition
				old_val = new_val
			end
		#else
		#	println("NOT VALID")
		end
		iter += 1
	end
	return partition, old_val
end
