using JuMP
using CPLEX
using CSV
using DataFrames

function dualisation(filename, sol_initiale=[])
	include(filename)
	l = zeros(Float64, n, n)
	for i in 1:n
		for j in 1:n
			l[i,j] = sqrt((coordinates[i,1] - coordinates[j,1])^2
							+(coordinates[i,2] - coordinates[j,2])^2)
		end
	end

	
	if sol_initiale != []
		startVal_x = zeros(n,n)
		startVal_y = zeros(n,K)
		for partie in sol_initiale
			for i in partie
				for j in partie
					startVal_x[i,j] = 1
				end
			end
		end
		for k in 1:size(sol_initiale)[1]
			for i in sol_initiale[k]
				startVal_y[i,k] = 1
			end
		end
	end
	
	# Create the model
	m = Model(CPLEX.Optimizer)
	set_optimizer_attribute(m, "CPXPARAM_TimeLimit", 1)
	
	### Variables
	# x[i, j] = 1 if (i, j) in same set
	@variable(m, x[i in 1:n, j in 1:n;i!=j], Bin)
	@variable(m, y[i in 1:n, k in 1:K], Bin)
	if sol_initiale != [] 
		println("build from heuristic")
		for i in 1:n
			for j in 1:n
				if i != j
					set_start_value(x[i,j], startVal_x[i,j])
				end
			end
			for k in 1:K 
				set_start_value(y[i,k], startVal_y[i,k])
			end
		end
	end
	@variable(m, mu[k in 1:K]>=0)
	@variable(m, lambda[i in 1:n, k in 1:K]>=0)
	@variable(m, alpha>=0)
	@variable(m, beta[i in 1:n, j in 1:n;i!=j]>=0)


	### Constraints
	#@constraint(m, y[1,1] == 1)

	
	@constraint(m, [k in 1:K], W*mu[k]
							   + sum(W_v[v]*lambda[v,k] for v in 1:n)
							   + sum(w_v[v]*y[v,k] for v in 1:n)
							   <= B)

	@constraint(m, [v in 1:n, k in 1:K], mu[k]+lambda[v,k]>=w_v[v]*y[v,k])
	@constraint(m, [i in 1:n, j in 1:n;i!=j], (alpha+beta[i,j])>=(lh[i]+lh[j])*x[i,j])
	
	@constraint(m, [i in 1:n], sum(y[i,k] for k in 1:K)==1)
	@constraint(m, [i in 1:n, j in 1:n, k in 1:K;i!=j], y[i,k]+y[j,k]-x[i,j]<=1)
	#@constraint(m, [i in 1:n, j in 1:n, k in 1:K;i!=j], y[i,k]-y[j,k]+x[i,j]<=1)
	#@constraint(m, [i in 1:n, j in 1:n, k in 1:K;i!=j], -y[i,k]+y[j,k]+x[i,j]<=1)
	
	### Objective
	@objective(m, Min, sum(l[i,j]*x[i,j] for i in 1:n, j in 1:n if i<j) + L*alpha+sum(3*beta[i,j] for i in 1:n for j in 1:n if i<j))


    ### Solve the problem
    start = time()
    optimize!(m)
    stop = time()
	filename_sol = "dualisation_sol.csv"
	if isfile(filename_sol)
		df = CSV.read(filename_sol, DataFrame)
	else
		df = DataFrame("nom_fichier"=> [], "time"=>[], "value"=>[], "sol"=>[])
	end

	no_sol = true
	sol = [[] for k in 1:K]
    for i in 1:n
        for k in 1:K
            if value(y[i,k])==1
                push!(sol[k], i)
				no_sol = false
            end
        end
    end

	print(filename, " ", array_to_string(sol))

	replace_row = false
	for row in eachrow(df)
		if row[:nom_fichier] == filename
			row[:time] = stop-start
			if no_sol
				row[:value] = "None"
				row[:sol] = "None"
			else
				row[:value] = JuMP.objective_value(m)
				row[:sol] = array_to_string(sol)
			end
			replace_row = true
		end
	end
	if !replace_row
		if no_sol
			push!(df, [filename, stop-start, "None", "None"])
		else
			push!(df, [filename, stop-start, JuMP.objective_value(m), array_to_string(sol)])
		end
	end

	CSV.write(filename_sol, df)
	return sol, JuMP.objective_value(m)
end


function array_to_string(array)
	print(array)
	result = ""
	for i in 1:size(array)[1]
		result = string(result, "[")
		for j in 1:length(array[i])
			if (j > 0)
				result = string(result, array[i][j])
			else
				result = string(result, array[i][j], ", ")
			end
		end
		result = string(result, "]")
	end
	return result
end