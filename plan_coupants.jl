function plan_coupants(filename, time_lim=60)
	include(filename)
	global l = zeros(Float64, n, n)
	for i in 1:n
		for j in 1:n
			l[i,j] = sqrt((coordinates[i,1] - coordinates[j,1])^2
							+(coordinates[i,2] - coordinates[j,2])^2)
	
		end
	end
	epsilon = 10^-6

	start = time()

	U1 = init_U1()
	U2 = init_U2()
	violation = true
	y_star,x_star, z_star, obj_lb, gap = -1,-1,-1,-1, -1
	while(violation)
		violation = false
		#println("ajout de contraintes")
		if time_lim - time()+start <1
			break
		end
		try
			x_star, y_star, z_star, obj_lb, gap = PM(U1, U2, time_lim - time()+start) 
		catch
			return
		end
		#println("y ", y_star)
		if time_lim - time()+start <1
			break
		end
		delta1_star = SP1(x_star, time_lim - time()+start)
		#println("z_star ",z_star)
		#println("sum ",sum(x_star[i,j]*(l[i,j]+delta1_star[i,j]*(lh[i]+lh[j])) for i in 1:n for j in 1:n if i<j))
		if(sum(x_star[i,j]*(l[i,j]+delta1_star[i,j]*(lh[i]+lh[j])) for i in 1:n for j in 1:n if i<j) - z_star > epsilon)
			#println("violation de contrainte sur les distances")
			l1 = zeros(Float64, n, n)
			for i in 1:n
				for j in 1:n
					l1[i,j] = l[i,j]+delta1_star[i,j]*(lh[i]+lh[j])
				end
			end            
			push!(U1,l1)
			violation = true
		end
		for k in 1:K
			if time_lim - time()+start <1
				break
			end
			delta2_star_k = SP2k(k, y_star, time_lim - time()+start)
			if(sum(y_star[i,k]*w_v[i]*(1+delta2_star_k[i]) for i in 1:n) - B > epsilon)
				#println("violation de contrainte sur les poids")
				w2 = Vector{Float64}(undef, n)
				for i in 1:n
					w2[i] = w_v[i]*(1+delta2_star_k[i])
				end
				push!(U2[k],w2)
				violation = true
			end
		end
	end
	stop = time()
	
	if violation != true
		sol = [[] for k in 1:K]
		for i in 1:n
			for k in 1:K
				if value(y_star[i,k])==1
					push!(sol[k], i)
				end
			end
		end
		write("plans_coupants", filename, stop - start, sol, z_star, obj_lb, gap)
		return sol, z_star
	end
	return
end

function init_U1()
	U1 = Vector{Matrix{Float64}}()
	push!(U1, l)
	return U1
end

function init_U2()
	U2 = Vector{Vector{Vector{Float64}}}(undef, K)
	for k in 1:K
		U2[k] = Vector{Vector{Float64}}()
		push!(U2[k], w_v)
	end
	return U2
end

function PM(U1, U2, time)
	# Create the model
	m = Model(CPLEX.Optimizer)
	set_optimizer_attribute(m, "CPX_PARAM_SCRIND", 0)
	set_silent(m)
	set_time_limit_sec(m, time)
	### Variables
	@variable(m, x[i in 1:n, j in 1:n; i!=j], Bin)
	@variable(m, y[i in 1:n, k in 1:K], Bin)
	@variable(m, z>=0)

	### Constraints
	@constraint(m, [u1 in U1], z>=sum(u1[i,j]*x[i,j] for i in 1:n for j in 1:n if i < j))
	@constraint(m, [k in 1:K, u2 in U2[k]], sum(u2[i]*y[i,k] for i in 1:n)<=B)
	
	@constraint(m, [i in 1:n], sum(y[i,k] for k in 1:K)==1)
	@constraint(m, [i in 1:n, j in 1:n, k in 1:K; i!=j], y[i,k]+y[j,k]-x[i,j]<=1)
	@constraint(m, [i in 1:n, j in 1:n, k in 1:K; i!=j], y[i,k]-y[j,k]+x[i,j]<=1)
	@constraint(m, [i in 1:n, j in 1:n, k in 1:K; i!=j], -y[i,k]+y[j,k]+x[i,j]<=1)
	
	### Objective
	@objective(m, Min, z)
	### Solve the problem
	optimize!(m)

	if(has_values(m))
		x_val = zeros(Int64, n, n)
		y_val = zeros(Int64, n, K)
		for i in 1:n
			for j in 1:n
				if i!=j
					x_val[i,j] = floor(value(x[i,j]))
				end
			end
			for k in 1:K
				y_val[i,k] = floor(value(y[i,k]))
			end
		end
		return x_val, y_val, value(z), objective_bound(m), relative_gap(m)
	else
		return
	end
end

function SP1(x_star, time)
	m = Model(CPLEX.Optimizer)
	set_optimizer_attribute(m, "CPX_PARAM_SCRIND", 0)
	set_silent(m)
	set_time_limit_sec(m, time)

	### Variables
	@variable(m, delta1[i in 1:n, j in 1:n; i!=j]>=0)

	### Constraints 
	@constraint(m, [i in 1:n, j in 1:n; i!=j], delta1[i,j]<=3)
	@constraint(m, sum(delta1[i,j] for i in 1:n, j in 1:n if i!=j)<=L)
	
	### Objective
	@objective(m, Max, sum((l[i,j]+(lh[i]+lh[j])*delta1[i,j])*x_star[i,j] for i in 1:n, j in 1:n if i!=j))
	
	### Solve the problem
	optimize!(m)
	delta_val = zeros(Float64, n, n)
	for i in 1:n
		for j in 1:n
			if i!=j
				delta_val[i,j] = value(delta1[i,j])
			end
		end
	end
	return delta_val
end

function SP2k(k, y_star, time)
	m = Model(CPLEX.Optimizer)
	set_optimizer_attribute(m, "CPX_PARAM_SCRIND", 0)
	set_silent(m)
	set_time_limit_sec(m, time)

	### Variables
	@variable(m, delta2[i in 1:n]>=0)

	### Constraints    
	@constraint(m, [i in 1:n], delta2[i] <= W_v[i])
	@constraint(m, sum(delta2[i] for i in 1:n)<=W)
	
	### Objective
	@objective(m, Max, sum(w_v[i]*y_star[i,k]*(1+delta2[i]) for i in 1:n))
	
	### Solve the problem
	optimize!(m)    
	return Vector{Float64}([value(delta2[i]) for i in 1:n])
end

