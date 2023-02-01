using JuMP
using CPLEX

function PL_statique(filename, time_lim)
	include(filename)
	l = zeros(Float64, n, n)
	for i in 1:n
		for j in 1:n
			l[i,j] = sqrt((coordinates[i,1] - coordinates[j,1])^2
							+(coordinates[i,2] - coordinates[j,2])^2)
	
		end
	end
	
	# Create the model
	m = Model(CPLEX.Optimizer)
	set_silent(m)
	set_time_limit_sec(m, time_lim)
	### Variables
	# x[i, j] = 1 if (i, j) in same set
		@variable(m, x[i in 1:n, j in 1:n; i!=j], Bin)
		@variable(m, y[i in 1:n, k in 1:K], Bin)

	### Constraints
	
	@constraint(m, [k in 1:K], sum(w_v[i]*y[i,k] for i in 1:n)<=B)
	
	@constraint(m, [i in 1:n], sum(y[i,k] for k in 1:K)==1)
	@constraint(m, [i in 1:n, j in 1:n, k in 1:K;i!=j], y[i,k]+y[j,k]<=x[i,j]+1)
	
	### Objective
	@objective(m, Min, sum(l[i,j]*x[i,j] for i in 1:n, j in 1:n if i<j))


	### Solve the problem
	start = time()
	optimize!(m)
	stop = time()
	sol = [[] for k in 1:K]
	for i in 1:n
		for k in 1:K
			if value(y[i,k])==1
				push!(sol[k], i)
			end
		end
	end
	write("statique", filename, stop-start, sol, objective_value(m), objective_bound(m), relative_gap(m))
end
