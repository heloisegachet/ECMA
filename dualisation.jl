using JuMP
using CPLEX
include("parser_out.jl")

function dualisation(filename, time_lim=60; gap=1e-6, sol_initiale=[], relax=nothing, cut_one=false)
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
	set_silent(m)
	set_time_limit_sec(m, time_lim)
	set_optimizer_attribute(m, "CPX_PARAM_EPGAP", gap)
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

	if !isnothing(relax)
		println("cut relax statique")
		@constraint(m, sum(l[i,j]*x[i,j] for i in 1:n, j in 1:n if i<j) + L*alpha+sum(3*beta[i,j] for i in 1:n for j in 1:n if i<j)>=relax)
	end
	if cut_one
		println("cut cluster one")
		@constraint(m,y[1,1]==1)
	end
	@constraint(m, [k in 1:K], W*mu[k]
							   + sum(W_v[v]*lambda[v,k] for v in 1:n)
							   + sum(w_v[v]*y[v,k] for v in 1:n)
							   <= B)

	@constraint(m, [v in 1:n, k in 1:K], mu[k]+lambda[v,k]>=w_v[v]*y[v,k])
	@constraint(m, [i in 1:n, j in 1:n;i!=j], (alpha+beta[i,j])>=(lh[i]+lh[j])*x[i,j])
	
	@constraint(m, [i in 1:n], sum(y[i,k] for k in 1:K)==1)
	@constraint(m, [i in 1:n, j in 1:n, k in 1:K;i!=j], y[i,k]+y[j,k]-x[i,j]<=1)
	@constraint(m, [i in 1:n, j in 1:n, k in 1:K;i!=j], y[i,k]-y[j,k]+x[i,j]<=1)
	@constraint(m, [i in 1:n, j in 1:n, k in 1:K;i!=j], -y[i,k]+y[j,k]+x[i,j]<=1)
	
	### Objective
	@objective(m, Min, sum(l[i,j]*x[i,j] for i in 1:n, j in 1:n if i<j) + L*alpha+sum(3*beta[i,j] for i in 1:n for j in 1:n if i<j))


    ### Solve the problem
    start = time()
    optimize!(m)
    stop = time()
	
	if(has_values(m))
		sol = [[] for k in 1:K]
		for i in 1:n
			for k in 1:K
				if value(y[i,k])==1
					push!(sol[k], i)
				end
			end
		end
		write("dualisation", filename, stop - start, sol, objective_value(m), string(objective_bound(m)), string(relative_gap(m)))
		
		return sol, objective_value(m), relative_gap(m)
	end
	return 
end

