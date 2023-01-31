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
	obj = 0
	for index_k in 1:size(partition)[1]
		if (length(partition[index_k]) >= 2)
			obj += sum(l_statique[i,j] + delta[i,j]*(lh[i]+lh[j]) for i in partition[index_k], j in partition[index_k] if i<j) 
		end
	end
	
	@objective(m, Max, obj)
	### Solve the problem
	start = time()
	optimize!(m)
	stop = time()
	#println(" solution of obj value ", JuMP.objective_value(m)," solving time = ",stop - start, "s")
	#for i in 1:n
	#	for j in 1:n
	#		println(i," ",j," l ", l_statique[i,j], " delta ", value(delta[i,j]), " lh ", lh[i]+lh[j])
	#	end
	#end
	return JuMP.objective_value(m)
end

function poids_partie(partie)
	# Create the model
	m = Model(CPLEX.Optimizer)
	set_silent(m)
	### Variables
	@variable(m, delta[v in partie]>=0)
	### Constraints
	@constraint(m, [v in partie], delta[v]<=W_v[v])
	@constraint(m, sum(delta[v] for v in partie) <=W)
	@objective(m, Max, sum(w_v[v]*(1+delta[v]) for v in partie))
	### Solve the problem
	start = time()
	optimize!(m)
	stop = time()
	somme = sum(w_v[v]*(1+value(delta[v])) for v in partie)
	#println("partie ", partie, " somme poids ", somme, " B ", B)
	#for v in partie
	#	println(" delta ", v, " ", value(delta[v]), " ", W_v[v], " ", W)
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
			break
		else
			# placer i dans la partie best_k et maj des valeurs
			push!(partition[best_k], i)
			assignment[i] = best_k
			w_group[best_k] = w_group[best_k] + w_nodes[i]
			v_group[best_k] += min_added_value
		end
	end

	println("RECHERCHE LOCALE")
	best_partition, best_val = recherche_locale(partition)

	stop = time()

	# ecrire la solution dans un fichier
	fout = open("./solution_heuristique.txt", "a")
	if v_sol == 0
		println(fout, "file "*filepath*" with n = ", n, "erreur", "solving time = ",stop - start, "s")
	else
		println(fout, "file "*filepath*" with n = ", n, " solution of obj value ", best_val, " solving time = ",stop - start, "s")
		println(fout, "solution : ")
		println(fout, best_partition)
	end
	close(fout)
	return best_partition, best_val
end


function recherche_locale(partition)
	old_val = value_sol(partition)
	iter = 0
	while iter < 100
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
			if poids_partie(partie) > B 
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
