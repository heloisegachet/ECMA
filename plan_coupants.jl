function plan_coupants()
    U1 = init_U1()
    U2 = init_U2()
    violation = true
    while(violiation)
        violation = false
        x_star, y_star, z_star = PM(U1, U2) 
        delta1_star = SP1(x_star)
        delta2_star = Vector{Float64}(0,3)
        delta2_star[k] = SP2k(k,y_star)
        if(sum(x_star[i,j]*(l[i,j]+delta1[i,j]*(lh[i]+lh[j]))) - z_star > epsilon)
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
            if(sum(y_star[i,k]*w_v[i]*(1+delta2_star[k][i]) - B > epsilon))
                w2 = Vector{Float64}(0, n)
                for i in 1:n
                    w2[i] = w_v[i]*(1+delta2_star[k][i])
                end
                push!(U2,w2)
                violation = true
            end
        end
end

function init_U1()
    U1 = Vector{Float64}()
    push!(U1, l)
    return U1
end

function init_U1()
    U2 = Vector{Float64}()
    push!(U2, w_v)
    return U2
end

function PM(U1, U2)
    # Create the model
    m = Model(CPLEX.Optimizer)

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
    start = time()
    optimize!(m)
end

function SP1(x_star)
    m1 = Model(CPLEX.Optimizer)

    ### Variables
    @variable(m1, 3>=delta1[i in 1:n, j in 1:n]>=0)

    ### Constraints    
    @constraint(m1, sum(delta1[i,j] for i in 1:n, j in 1:n)<=L)
    
    ### Objective
    @objective(m1, Max, sum(l[i,j]*x_star[i,j]+(lh[i]+lh[j])*delta1[i,j]*x_star[i,j] for i in 1:n, j in 1:n))
    
    ### Solve the problem
    optimize!(m1)
    return JuMP.value(delta1)
end

function SP2k(k, y_star)
    m2 = Model(CPLEX.Optimizer)

    ### Variables
    @variable(m2, W_v[i]>=delta2[i in 1:n]>=0)

    ### Constraints    
    @constraint(m2, sum(delta2[i] for i in 1:n)<=W)
    
    ### Objective
    @objective(m2, Max, sum(w_v[i]*y_star[i,k]+w_v[i]*delta2[i]*y_star[i,k] for i in 1:n))
    
    ### Solve the problem
    optimize!(m2)
    return JuMP.value(delta2)
end