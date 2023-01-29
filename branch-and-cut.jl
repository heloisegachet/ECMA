using JuMP
using CPLEX


function isIntegerPoint(cb_data::CPLEX.CallbackContext, context_id::Clong)
	# context_id == CPX_CALLBACKCONTEXT_CANDIDATE si le callback est
	# appelé dans un des deux cas suivants :
	# cas 1 - une solution entière a été obtenue; ou
	# cas 2 - une relaxation non bornée a été obtenue
	if context_id != CPX_CALLBACKCONTEXT_CANDIDATE
		return false
	end
	# Pour déterminer si on est dans le cas 1 ou 2, on essaie de récupérer la
	# solution entière courante
	ispoint_p = Ref{Cint}()
	ret = CPXcallbackcandidateispoint(cb_data, ispoint_p)
	# S’il n’y a pas de solution entière
	if ret != 0 || ispoint_p[] == 0
		return false
	else
		return true
	end
end

function branch_and_cut(filepath)
	include(filepath)
	global l = zeros(Float64, n, n)
	for i in 1:n
		for j in 1:n
			l[i,j] = sqrt((coordinates[i,1] - coordinates[j,1])^2
							+(coordinates[i,2] - coordinates[j,2])^2)
		end
	end
	epsilon = 10^-6

	# Create the model
	m = Model(CPLEX.Optimizer)
	set_optimizer_attribute(m, "CPX_PARAM_SCRIND", 0)
	MOI.set(m, MOI.NumberOfThreads(), 1)

	### Variables
	@variable(m, x[i in 1:n, j in 1:n; i!=j], Bin)
	@variable(m, y[i in 1:n, k in 1:K], Bin)
	@variable(m, z>=0)

	### Constraints	
	@constraint(m, z>=sum(l[i,j]*x[i,j] for i in 1:n for j in 1:n if i < j))
	@constraint(m, [k in 1:K], sum(w_v[i]*y[i,k] for i in 1:n)<=B)
	
	@constraint(m, [i in 1:n], sum(y[i,k] for k in 1:K)==1)
	@constraint(m, [i in 1:n, j in 1:n, k in 1:K; i!=j], y[i,k]+y[j,k]-x[i,j]<=1)
	@constraint(m, [i in 1:n, j in 1:n, k in 1:K; i!=j], y[i,k]-y[j,k]+x[i,j]<=1)
	@constraint(m, [i in 1:n, j in 1:n, k in 1:K; i!=j], -y[i,k]+y[j,k]+x[i,j]<=1)

	### Objective
	@objective(m, Min, z)

	### Callback function
	function lazy_callback(cb_data::CPLEX.CallbackContext, context_id::Clong)
		if isIntegerPoint(cb_data, context_id)
			# Cette ligne doit être appelée avant de pouvoir récupérer la
			# solution entière ayant entraîné l’appel du callback
			CPLEX.load_callback_variable_primal(cb_data, context_id)
			# On récupère la valeur de x, y, z
			z_star = callback_value(cb_data, z)
			x_star = zeros(Int64, n, n)
			#print(z_star)
			for i in 1:n
				for j in 1:n
					if j != i
						x_star[i,j] = round(callback_value(cb_data, x[i,j]))#round(callback_value(cb_data, m[:x][i,j]))
					end
				end
			end
			y_star = zeros(Int64, n, K)
			for i in 1:n
				for k in 1:K
					y_star[i,k] = round(callback_value(cb_data, y[i,k]))#round(callback_value(cb_data, m[:y][i,k]))
				end
			end
			
			### Find a U1 cut
			delta1_star = SP1(x_star)
			if(sum(x_star[i,j]*(l[i,j]+delta1_star[i,j]*(lh[i]+lh[j])) for i in 1:n for j in 1:n if i<j) - z_star > epsilon)
				#println("cut the length")
				cstr1 = @build_constraint(z>=sum(x[i,j]*(l[i,j]+delta1_star[i,j]*(lh[i]+lh[j])) for i in 1:n for j in 1:n if i<j))
				MOI.submit(m, MOI.LazyConstraint(cb_data), cstr1)
			end
			### Find a U2 cut
			for k in 1:K
				delta2_star_k = SP2k(k, y_star)
				if(sum(y_star[i,k]*w_v[i]*(1+delta2_star_k[i]) for i in 1:n) - B > epsilon)
					#println("cut the weight")
					cstr2 = @build_constraint(sum(y[i,k]*w_v[i]*(1+delta2_star_k[i]) for i in 1:n)<=B)
					MOI.submit(m, MOI.LazyConstraint(cb_data), cstr2)
				end
			end
		end
	end

	### Using callback function
	MOI.set(m, CPLEX.CallbackFunction(), lazy_callback)

	### Solve
	start = time()
	optimize!(m)
	stop = time()

	### Check
	println(JuMP.objective_value(m))
	println(stop-start, "s")
end


function SP1(x_star)
	m = Model(CPLEX.Optimizer)
	set_optimizer_attribute(m, "CPX_PARAM_SCRIND", 0)

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
				delta_val[i,j] = JuMP.value(delta1[i,j])
			end
		end
	end
	return delta_val
end


function SP2k(k, y_star)
	m2 = Model(CPLEX.Optimizer)
	set_optimizer_attribute(m2, "CPX_PARAM_SCRIND", 0)

	### Variables
	@variable(m2, delta2[i in 1:n]>=0)

	### Constraints    
	@constraint(m2, [i in 1:n], delta2[i] <= W_v[i])
	@constraint(m2, sum(delta2[i] for i in 1:n)<=W)
	
	### Objective
	@objective(m2, Max, sum(w_v[i]*y_star[i,k]*(1+delta2[i]) for i in 1:n))
	
	### Solve the problem
	optimize!(m2)    
	return Vector{Float64}([JuMP.value(delta2[i]) for i in 1:n])
end