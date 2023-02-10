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
	@objective(m, Max, sum(sum(l[i,j] + delta[i,j]*(lh[i]+lh[j]) for i in partition[index_k], j in partition[index_k] if i<j) for index_k in 1:size(partition)[1] if (length(partition[index_k]) >= 2)))
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

function heuristic(filename, time_lim; val_min=nothing)

	include(filename)

	start = time()

	# liste des pires poids des sommets
	w_nodes = Vector{Float64}(zeros(n))
	for i in 1:n
		w_nodes[i] = w_v[i]*(1 + W_v[i])
	end

	# matrice des pires longueurs des arêtes
	global l = zeros(Float64, n, n)
	l_robuste = zeros(Float64, n, n)

	for i in 1:n
		for j in 1:n
			if j!=i
				l[i,j] = sqrt((coordinates[i,1] - coordinates[j,1])^2 +(coordinates[i,2] - coordinates[j,2])^2)
				l_robuste[i,j] = l[i,j] +3*(lh[i]+lh[j])
			end
		end
	end

	# constante de majoration (somme des longueurs)
	M = 0
	for i in 1:n
		for j in i+1:n
			M = M + l_robuste[i,j]
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
				v = added_value(i, partition[k], l_robuste)
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
	best_partition, best_val, best_time, iter = recherche_locale(partition, time_lim - time() + start, val_min)

	stop = time()

	write("heuristique", filename,best_time-start, best_partition, string(old_val," / ",best_val), string("iter : ",iter), "None")
	return best_partition, best_val
end


function recherche_locale(partition, time_lim, val_min)
	start = time()
	old_val = value_sol(partition)
	best_time = time()
	iter = 0
	iter = 0
	epsilon = 1e-10
	println(val_min, " ", old_val)
	while time() - start < time_lim
		if !isnothing(val_min) && old_val <= val_min + epsilon
			println(old_val, " , ", val_min + epsilon)
			println("optimum")
			break
		end
		iter += 1
		k1 = rand(1:K)
		k2 = rand([1:k1-1;k1+1:K])
		v1 = rand([partition[k1];-1])
		v2 = rand(partition[k2])
		partie1 = Vector{Int64}()
		partie2 = Vector{Int64}()
		for i in partition[k1]
			if i != v1
				push!(partie1, i)
			end
		end
		for i in partition[k2]
			if i != v2
				push!(partie2, i)
			end
		end
		if v1 != -1
			push!(partie2, v1)
		end
		push!(partie1, v2)

		if length(partie2)>=1 && sum(w_v[v] for v in partie2) > B
			#println("	first verif poids")
			continue
		end	
		if length(partie1)>=1 && sum(w_v[v] for v in partie1) > B
			#println("	first verif poids")
			continue
		end	

		if obj_statique(partition, partie1, k1, partie2, k2) > old_val
			#println("	first verif obj")
			continue
		end
		


		poids1 = poids_partie(partie1, time_lim - time()+start)
		if isnothing(poids1)
			println("iter ",iter)
			return partition, old_val, best_time, iter
		end
		if poids1 > B 
			#println("	verif poids PL")
			continue
		end

		poids2 = poids_partie(partie2, time_lim - time()+start)
		if isnothing(poids2)
			println("iter ",iter)
			return partition, old_val, best_time, iter
		end
		if poids2 > B 
			#println("	verif poids PL")
			continue
		end

		#println("		test valeur echange ")
		new_val = value_before_change(partition, k1, partie1, k2, partie2, time_lim-time()+start)
		if old_val > new_val
			#println("	better replace : ", new_val)
			partition[k1] = copy(partie1)
			partition[k2] = copy(partie2)
			old_val = new_val
			best_time = time()
		end
		#else
		#	println("NOT VALID")
	end
	println("iter ",iter)
	return partition, old_val, best_time, iter
end


function value_before_change(partition, k1, partie1, k2, partie2, time_lim)
	# Create the model
	m = Model(CPLEX.Optimizer)
	set_silent(m)
	set_time_limit_sec(m, time_lim)

	### Variables
	@variable(m, delta[i in 1:n, j in 1:n]>=0)
	### Constraints
	@constraint(m, [i in 1:n, j in 1:n;i!=j], delta[i,j]<=3)
	@constraint(m, sum(delta[i,j] for i in 1:n, j in 1:n if i!=j) <=L)
	### Objective
	
	@objective(m, Max, sum(
							sum(l[i,j] + delta[i,j]*(lh[i]+lh[j]) for i in partition[index_k], j in partition[index_k] if i<j) 
							for index_k in 1:size(partition)[1] if (length(partition[index_k]) >= 2 && index_k != k1 && index_k != k2)
							)
						+ sum(l[i,j] + delta[i,j]*(lh[i]+lh[j]) for i in partie1, j in partie1 if i<j)
						+ sum(l[i,j] + delta[i,j]*(lh[i]+lh[j]) for i in partie2, j in partie2 if i<j)
			    )
	### Solve the problem
	optimize!(m)

	#for index_k in 1:size(partition)[1] 
	#	println("partie ", partition[index_k])
	#	if (length(partition[index_k]) >= 2 && index_k != k1 && index_k != k2)
	#		println("val ",sum(l[i,j] + value(delta[i,j])*(lh[i]+lh[j]) for i in partition[index_k], j in partition[index_k] if i<j))
	#	end
	#end
	#println("partie ", partition[k1])
	#if (length(partition[k1]) >= 2)
	#	println("val ",sum(l[i,j] + value(delta[i,j])*(lh[i]+lh[j]) for i in partition[k1], j in partition[k1] if i<j))
	#end
	#println("partie ", partition[k2])
	#if (length(partition[k2]) >= 2)
	#	println("val ",sum(l[i,j] + value(delta[i,j])*(lh[i]+lh[j]) for i in partition[k2], j in partition[k2] if i<j))
	#end

	return objective_value(m)
end

function obj_statique(partition, partie1, k1, partie2, k2)
	somme = 0
	if length(partie1)>=2
		somme += sum(l[i,j] for i in partie1 for j in partie1 if i<j)
	end
	if length(partie2)>=2
		somme += sum(l[i,j] for i in partie2 for j in partie2 if i<j)
	end
	for k in size(partition)[1]
		if k != k1 && k!= k2
			if length(partition[k]) >= 2
				somme += sum(l[i,j] for i in partition[k] for j in partition[k] if i<j)
			end
		end
	end
	return somme
end