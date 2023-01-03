function plan_coupants(filename)
    include(filename)
    l = zeros(Float64, n, n)
    for i in 1:n
        for j in 1:n
            l[i,j] = sqrt((coordinates[i,1] - coordinates[j,1])^2
                            +(coordinates[i,2] - coordinates[j,2])^2)
    
        end
    end
    epsilon = 10^-6

    U1 = init_U1(l)
    U2 = init_U2()
    violation = true
    while(violation)
        violation = false
        x_star, y_star, z_star = PM(U1, U2,l) 
        delta1_star = SP1(x_star, l)
        
        if(sum(x_star[i,j]*(l[i,j]+delta1_star[i,j]*(lh[i]+lh[j])) for i in 1:n for j in 1:n) - z_star > epsilon)
            println("violation de contrainte sur les distances")
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
            delta2_star = SP2k(k, y_star)
            if(sum(y_star[i,k]*w_v[i]*(1+delta2_star[i]) for i in 1:n) - B > epsilon)
                println("violation de contrainte sur les poids")
                w2 = Vector{Float64}(undef, n)
                for i in 1:n
                    w2[i] = w_v[i]*(1+delta2_star[i])
                end
                push!(U2,w2)
                violation = true
            end
        end
    end
end

function init_U1(l)
    U1 = Vector{Matrix{Float64}}()
    push!(U1, l)
    return U1
end

function init_U2()
    U2 = Vector{Vector{Float64}}()
    push!(U2, w_v)
    return U2
end

function PM(U1, U2, l)
    # Create the model
    m = Model(CPLEX.Optimizer)
    set_optimizer_attribute(m, "CPX_PARAM_SCRIND", 0)

    ### Variables
    @variable(m, x[i in 1:n, j in 1:n], Bin)
    @variable(m, y[i in 1:n, k in 1:K], Bin)
    @variable(m, z>=0)

    ### Constraints
    @constraint(m, [l1 in U1], z>=sum(l1[i,j]*x[i,j] for i in 1:n for j in 1:n))
    @constraint(m, [k in 1:K, w2 in U2], sum(w2[i]*y[i,k] for i in 1:n)<=B)
    
    @constraint(m, [i in 1:n], sum(y[i,k] for k in 1:K)==1)
    @constraint(m, [i in 1:n, j in 1:n, k in 1:K], y[i,k]+y[j,k]<=x[i,j]+1)
    
    ### Objective
    @objective(m, Min, sum(l[i,j]*x[i,j] for i in 1:n, j in 1:n if i<j))
    ### Solve the problem
    optimize!(m)

    x_val = zeros(Int64, n, n)
    y_val = zeros(Int64, n, K)
    for i in 1:n
        for j in 1:n
            x_val[i,j] = JuMP.value(x[i,j])
        end
        for k in 1:K
            y_val[i,k] = JuMP.value(y[i,k])
        end
    end

    return x_val, y_val, JuMP.value(z)
end

function SP1(x_star, l)
    m1 = Model(CPLEX.Optimizer)
    set_optimizer_attribute(m1, "CPX_PARAM_SCRIND", 0)

    ### Variables
    @variable(m1, 3>=delta1[i in 1:n, j in 1:n]>=0)

    ### Constraints    
    @constraint(m1, sum(delta1[i,j] for i in 1:n, j in 1:n)<=L)
    
    ### Objective
    @objective(m1, Max, sum(l[i,j]*x_star[i,j]+(lh[i]+lh[j])*delta1[i,j]*x_star[i,j] for i in 1:n, j in 1:n))
    
    ### Solve the problem
    optimize!(m1)
    delta_val = zeros(Float64, n, n)
    for i in 1:n
        for j in 1:n
            delta_val[i,j] = JuMP.value(delta1[i,j])
        end
    end
    return delta_val
end

function SP2k(k, y_star)
    m2 = Model(CPLEX.Optimizer)
    set_optimizer_attribute(m2, "CPX_PARAM_SCRIND", 0)

    ### Variables
    @variable(m2, W_v[i]>=delta2[i in 1:n]>=0)

    ### Constraints    
    @constraint(m2, sum(delta2[i] for i in 1:n)<=W)
    
    ### Objective
    @objective(m2, Max, sum(w_v[i]*y_star[i,k]+w_v[i]*delta2[i]*y_star[i,k] for i in 1:n))
    
    ### Solve the problem
    optimize!(m2)    
    return Vector{Float64}([JuMP.value(delta2[i]) for i in 1:n])
end