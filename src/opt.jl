#using LinearAlgebra, Plots, ForwardDiff, DataFrames
#pyplot(framestyle = :box, color = :black, linestyle = :auto)

# Assignment 1

## Descent Functions for gradient-based optimization

### Steepest Descent
function descent_steepest(f::Function, position_current::Vector)
    # Determine the gradient at selected position
    grad_current = ForwardDiff.gradient(f, position_current)
    
    # Steepest descent
    d = - grad_current
    
    # Output descent direction (vector), gradient at current step (Vector)
    return d
end

### Conjugate gradient
function descent_conjugate(grad_current::Vector, grad_prev::Vector, d_prev::Vector)
    #Determine the value of β_k for the Conjugate Gradient search method
    β = (norm(grad_current) / norm(grad_prev))^2
    
    #Determine new direction vector for current step
    d = -grad_current .+ β * d_prev
    
    # Output descent direction (Vector), gradient at current step (Vector)
    return d
end

### Newton
function descent_newton(grad_current::Vector, hess_current::Array, position_current::Vector; det_err = 1e-8)
    #Two conditions must be met for newton method to be valid.
    
    #Condition 1: Positive Semi-Definite Hessian
    #IE if there is an eigenvalue of the Hessian that is less than 0, the Hessian is not positive semi-definite
    #and the newton method cannot be used
    if any(eigvals(hess_current) .< 0)
        println("Hessian of function at this point is not positive semi-definite. Newton method cannot be used.")
        return
    end
    
    #Condition 2: Hessian must be non-singular
    #IE: The Hessian must be invertible,
    #IE2: The determinant of the Hessian must be non-zero
    #Note that det_err provides a safeguard against numerical errors when checking for the determinant value
    if abs(det(hess_current)) < det_err
        println("Hessian of function at this point is singular. Newton method cannot be used.")
        return
    end
    
    #If both requirements pass, then the newton method provides the output:
    d = -inv(hess_current) * grad_current
    
    return d
end

# Here's another git test 

## Step Function (Line search)
function α_step(func::Function, position_current::Vector, d_current::Vector, α_bar; ϵ = 1e-4, intervals = 50, max_iter = 100)
    #α_bar is the upper range for α, set as 1.0 typically, IE 100% of the gradient
    
    #Initiate values
    err = 1 
    lower = 0 #Lower bound is 0, giving f(x_(k+1)) = f(x_k)
    upper = α_bar #Upper bound is the user defined α_bar
    iter = 0
    
    while err > ϵ #While the error remains large
        
        #Update iteration count
        iter += 1
        
        if iter > max_iter
            println("Iteration time out for search. Increase maximum iteration number or revise input.")
            return
        end
        
        #Range of increasing values of α from 0 to α_bar
        α_range = collect(range(lower, length = intervals, stop = upper))
        #This creates an array of length = 'intervals' from lower bound to upper bound
    
        #f_test_range is the value of f(x) for each increment α_i ∈ α_range
        f_test_range = [func(position_current .+ α_range[i] .* d_current) for i = 1:intervals]
        
        #Find the incremental difference in value of f_test_range
        f_diff = diff(f_test_range)
        
        #Note that it is possible that the selected range has a continuously decreasing function, in which case
        #the optimum value of α is not determined in the selection range, therefore incerase the upper bound on α
        if all(f_diff .< 0)
            upper = 2 * upper
            continue
        end
        
        #find the index of the first instance where the function begins to increase
        idx_lower = findfirst(x-> x>0, f_diff) - 1 #lower bound of new interval
        idx_upper = idx_lower + 2 #upper bound of new interval
    
        #Updated lower and upper bounds for the optimal step size α
        lower = α_range[idx_lower]
        upper = α_range[idx_upper]
        
        #Calculate error
        err = upper - lower
        
    end
    
    #After this loop, the interval for the upper and lower limits for α has been determined
    #Return the average between the bounds
    α_opt = (upper + lower) / 2
    
    #return the optimal step size + number of iterations taken to reach it
    return α_opt, iter

end

## Main Unconstrainted Optimization function
#Define the error calculating function
function gradnorm(func::Function, x::Vector)
    return norm(ForwardDiff.gradient(func, x))
end

function opt_unconstrained(func::Function, 
        x_start::Array,
        descent_method,
        iter_max::Int, 
        α_max::Float64,
        ϵ::Float64;
        annotation = true)
    
    #Initiate storage vectors for relevant data:
    iter_store = [] #number of iterations 
    x_store = [] #value of design variables at each iteration
    f_store = [] #value of function at each iteration
    grad_store = [] #gradient at each iteration
    err_store = [] #error value at each iteration (norm of gradient)
    iter_linesearch_store = [] #number of iterations taken for line search algorithm
    α_store = [] #step size
    d_store = [] #descent vector
    
    #initiate values value
    err = 1
    x = x_start
    iter = 0
    
    while err > ϵ && iter < iter_max
        iter += 1
        
        #Iteration values
        f = func(x)
        grad = ForwardDiff.gradient(func,x)
        hess = ForwardDiff.hessian(func,x)
        err = gradnorm(func, x)
        
        #Find descent direciton, d
        
        if descent_method == :steepest
            d = -grad
        #NOTE:Since the conjugate gradient method relies on the gradient and direction vector
        # of the previous step, the direction vector for the FIRST step is set to the steepest descent method
        elseif descent_method == :cg && iter == 1
            d = -grad
        elseif descent_method == :cg
            d = descent_conjugate(grad, grad_store[iter-1], d_store[iter-1])
        elseif descent_method == :newton
            d = descent_newton(grad, hess, x)
        else
            println("Descent method unrecognized.")
            println("Allowable methods: :steepest, :cg, :newton")
            return
        end
            
        #find the step size, α, and the number of iterations to solve for it, α_iter
        α, α_iter = α_step(func, x, d, α_max)
        
        #Store all relevant values
        push!(iter_store, iter -1) #number of main loop iterations
        push!(x_store, x) #value of design variables
        push!(f_store, f) #value of function
        push!(grad_store, grad) #value of gradient
        push!(err_store, err) #value of error (norm of gradient)
        push!(α_store, α) #value of step size
        push!(iter_linesearch_store, α_iter) #number of line search iterations
        push!(d_store, d) #step direction vector
        
        #update the value of the design variables, 
        x = x .+ α .* d
        
        if annotation
            println("Iteration: ", iter-1)
            println("Function value: ", round(f, digits = 2))
            println("Error (Norm(Gradient)): ", round(err, digits = 2))
            println("Line search iterations: ", α_iter)
            println("Step size: ", round(α, digits = 2))
            println("________________________")
        end
    end

    if iter == iter_max
        println("Iteration time out. Revise.")
        return
    end
    
    #Output all necessary information as a DataFrame type for easy export + plotting
    output = DataFrame(N_iteration = iter_store,
        x_val = x_store,
        Function_val = f_store,
        Gradient = grad_store,
        Gradient_norm = err_store,
        Direction_vec = d_store,
        Step_size = α_store,
        Linesearch_iterations = iter_linesearch_store)
    
    return output
end

# Gradient computation

function adjoint_gradient(nodes, elements, loads)
    
    ## Initialize basic information
    lengths, angles, dof_active, dof_list, n_dof = fem_init(nodes, elements)
    
    #initialize force vector
    Force_vector = loadvector(loads, nodes, n_dof)
    
    #Create local/global elemental stiffness matrices
    T = [Γ(angle) for angle in angles]
    k_element_local = [k_localxy(elements[i], lengths[i]) for i = 1:length(elements)]
    k_element_global = [k_globalxy(k_element_local[i], T[i]) for i = 1:length(elements)]
    
    #Create global stiffness matrix
    K_global = build_K(n_dof, dof_active, nodes, elements, k_element_global)
    
    ##############################
    #Stiffness partial derivative matrix
    #############################
    K_sensitivity = []
    for i = 1:length(elements)
        K_temp = zeros(n_dof, n_dof)
        #for element i, determine which nodes are connected to i, and at what end
        n_start, n_end = start_end_node_idx(elements[i], nodes)
        
        #Find the index values for the active DOF in the global XY element stiffness matrix
        idx_local_start = [1,2] .* dof_active[n_start]
        idx_local_start = idx_local_start[idx_local_start .!= 0]
    
        #Find the index values for the active DOf in the global XY element stiffness matrix
        idx_local_end = [3,4] .* dof_active[n_end]
        idx_local_end = idx_local_end[idx_local_end .!= 0]
        
        #Concatenate the activated index values for the global XY element stiffness matrix
        idx_local = vcat(idx_local_start, idx_local_end)
        
        #Convert the above calculated values to the index of the GLOBAL stiffness matrix
        idx_global_start =  (n_start * 2 - 2) .+ idx_local_start
        idx_global_end =  (n_end * 2 - 4) .+ idx_local_end
        idx_global = vcat(idx_global_start, idx_global_end)
        
        #Extracted sub matrix from the global XY element stiffness matrix
        k_extracted = k_element_global[i][idx_local, idx_local] ./ elements[i].A
        
        #This is the global (Full DOF, active/inactive) stiffness matrix w/r/t element i
        K_temp[idx_global, idx_global] .+= k_extracted
        
        push!(K_sensitivity, K_temp)

    end 
    ##############################
    #Determine Compliance
    #############################
    
    F_reduced = Force_vector[dof_list]
    K_glob_reduced = K_global[dof_list, dof_list]
    
    disp = K_glob_reduced \ F_reduced
    
    K_sens_reduced = [K_sen[dof_list, dof_list] for K_sen in K_sensitivity]
    
    grad_adj = [- transpose(disp) * K_sen * disp for K_sen in K_sens_reduced]
    return grad_adj
end

function fd_gradiant(func, x_init, step)
    n_dof = length(x_init)
    f_init = func(x_init)
    grad = []
    for i = 1:n_dof
        x_step = zeros(n_dof)
        x_step[i] = step
        
        f_step = func(x_init + x_step)
        grad_step = (f_step - f_init) / step
        push!(grad, grad_step)
    end
    return grad
end
        
