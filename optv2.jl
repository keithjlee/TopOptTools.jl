using LinearAlgebra, Plots, NLopt, SparseArrays

# Geometric Functions
#Get length of element
function elem_length(element)
    return norm(element[2] - element[1])
end

function element_maker(nodes, element_indices)
    return [[nodes[e_idx[1]], nodes[e_idx[2]]] for e_idx in element_indices]
end


# get angle of element
function θ(element)
    x_global = [1, 0] #unit x vector
    vec_element = element[2] - element[1]

    cos = dot(x_global, vec_element) / (norm(x_global) * norm(vec_element))

    angle = acosd(cos)

    if vec_element[2] < 0
        angle = 360 - angle
    end

    return angle
end

## Transformation matrix
function Γ(ang)
    return [cosd(ang) sind(ang) 0 0;0 0 cosd(ang) sind(ang)]
end

function e_dof_expanded(e_idx)
    x_idx = e_idx .* 2 .- 1
    y_idx = x_idx .+ 1
    
    edof = zeros(Int, 4)
    edof[1:2:end-1] = x_idx
    edof[2:2:end] = y_idx
    
    return edof
end


## Loads
function loadidx(pos, nodes)
    # will create a load value to the closest node at the chosen position
    # pos = [x, y]
    
    load_idx = findmin([sum(abs.(n .- pos)) for n in nodes])[2]
    
    return load_idx
end

function loadmaker(loads, positions, nodes, n_dofs)
    if length(loads) !== length(positions)
        println("load and position vectors must be of equal length.")
        return
    end
    #find the nodal index for each load
    idx = [loadidx(pos, nodes) for pos in positions]
    idx_expanded = vcat([[i*2-1, i*2] for i in idx]...)
    loads_expanded = vcat(loads...)
    
    f = zeros(n_dofs)
    f[idx_expanded] = loads_expanded
    return f
end

##Stiffness

## Local coordinate elemental stiffness
function k_localxy(length, A, E)
    k = A * E / length .* [1 -1; -1 1]
    return k
end

# Convert local stiffness matrices to global stiffness matrices
function k_globalxy(k_local::Array{Float64, 2}, T)
    #Global_xy stiffness matrix
    k = transpose(T) * k_local * T
    
    #Note the indices of the matrix are in the form of:
    # 1: node A, global X
    # 2: node A, global Y
    # 3: node B, global X
    # 4: node B, global Y

    return k
end

# Assemble global stiffness matrix
function build_K(element_node_idx, k_elemental, n)
    # Indices in global stiffness matrix
    idxx = vcat([vcat([ei * ones(Int, 4) for ei in eidx]...) for eidx in element_node_idx]...) 
    idxy = vcat([repeat(ei, 4) for ei in element_node_idx]...)

    kvals = vcat([vcat(reshape(k_elemental[i], (16,1))...) for i = 1:n]...)

    return sparse(idxx, idxy, kvals)
end

## Solving

function displacements(K, F, dof_list, n_dof)
    #Reduce K matrix and force vector to non-zero items:
    K_reduced = K[dof_list, dof_list]
    F_reduced = F[dof_list]
    
    disp = Symmetric(K_reduced) \ F_reduced


    compliance = transpose(F_reduced) * disp #F^T d 
    
    #Extended disp vector
    disp_long = zeros(n_dof)
    disp_long[dof_list] = disp
    
    return disp_long, compliance
end

##Post Processing
function reactions(F_e, n_dofs, element_idx)
    #create empty storage matrix
    F_all = zeros(length(dofs))
    
        #for each element
        for i = 1:length(elements)
            #extract the global XY elemental end forecs
            F = F_e[i]
    
            #extract the given element information
            e = element_idx[i]
    
            #Determine the start and end nodes for the given element
            n_start = e[1]
            n_end = e[2]
    
            #Convert to the indexes in the global DOF list 
            idx_start = 2 * n_start - 1 .+ [0,1]
            idx_end = 2 * n_end - 1 .+ [0,1]
    
            #Add the end forces of each element to the appropriate DOF:
            F_all[idx_start] .+= F[1:2]
            F_all[idx_end] .+= F[3:4]
        end
        
    #Return the force values of each DOF that are at REACTIONS
    return F_all .* (dofs .== 0)
end

function forces(n_elements, elements_idx, A, k_global, disp, T, n_dofs)
    #Global element forces
    F_e = [k_global[i] * disp[i] for i = 1:n_elements]
    
    #Convert global element forces t olocal element forces
    f_e = [T[i] * F_e[i] for i = 1:n_elements]
    
    #Single axial force value
    f = [force[2] for force in f_e]
    
    #Stress values
    σ = [f[i] / A[i] for i = 1:n_elements]
    
    #Determing Reaction Forces
    rxns = reactions(F_e, n_dofs, elements_idx)
    
    return F_e, f_e, f, σ, rxns
end

function eq_check(rxns, F, tolerance)
    #Check equilibrium
    if sum(rxns) + sum(F) > tolerance
        println("Equilibrium not met!!")
        return true
    else
        println("Equilibrium met.")
    end
end

function eq_check(rxns, F, tolerance)
    #Check equilibrium
    if sum(rxns) + sum(F) > tolerance
        println("Equilibrium not met!!")
        return true
    else
        println("Equilibrium met.")
    end
end


function displaced_nodes(nodes, disp; scale_x = 1, scale_y = 1)
    n_disp = []
    for i = 1:length(nodes)
        n = nodes[i]
        n_new = [(n[1] + disp[i*2 - 1]) * scale_x, (n[2] + disp[i*2]) * scale_y]
        push!(n_disp, n_new)
    end
    return n_disp
end

## Secondary
function fem_init(nodes, elements)
    #extract lengths of each element
    lengths = [elem_length(e) for e in elements]
    #extract CC angle (degrees) of each element from global axis
    T = [Γ(θ(e)) for e in elements]
    
    return lengths, T
end

function stiffness_init(nodes, elements, lengths, T, A, E)

    #Stiffness Matrices
    k_element_local = [k_localxy(lengths[i], A[i], E[i]) for i = 1:length(elements)]
    k_element_global = [k_globalxy(k_element_local[i], T[i]) for i = 1:length(elements)]
    
    return k_element_global
end

## Main

## Main solver
function fem2d_solver(nodes, dofs, element_idx, A, E, loads, positions; tol = 1e-5)
    
    elements = element_maker(nodes, element_idx)
    n_elements = length(elements)
    lengths = [elem_length(e) for e in elements]
    
    T = [Γ(θ(e)) for e in elements]
    
    n_dof = length(dofs)

    #Return global stiffness matrix (including 0 rows/columns)
    k_global = stiffness_init(nodes, elements, lengths, T, A, E)
    
    idx_expanded = e_dof_expanded.(eidx)
    
    K = build_K(idx_expanded, k_global, n_elements)
    #Create force vector
    F = loadmaker(loads, positions, nodes, n_dof)
    
    #Displacements in global coordinate system + elemental DOF displacements
    disp, compliance = displacements(K, F, dofs, n_dof)
    
    disp_local = [disp[idx] for idx in idx_expanded]
    
    F_global, f_elemental, f_axial, stress_axial, rxns = forces(n_elements, element_idx, A, k_global, disp_local, T, n_dof)
    
    #If equilibrium is not reached, exit solver
    if eq_check(rxns, F, tol) == true
        return
    end
    
    #Deformed nodes and elements
    new_nodes = displaced_nodes(nodes, disp)
    
    return new_nodes, disp, f_axial, stress_axial, rxns, compliance
end

## Main solver
function analysis(nodes, dofs, element_idx, A, E, loads, positions; tol = 1e-5)
    
    elements = element_maker(nodes, element_idx)
    n_elements = length(elements)
    lengths = [elem_length(e) for e in elements]
    
    T = [Γ(θ(e)) for e in elements]
    
    n_dof = length(dofs)

    #Return global stiffness matrix (including 0 rows/columns)
    k_global = stiffness_init(nodes, elements, lengths, T, A, E)
    
    idx_expanded = e_dof_expanded.(eidx)
    
    K = build_K(idx_expanded, k_global, n_elements)

    if det(K) < 0
        return 0
    end
    #Create force vector
    F = loadmaker(loads, positions, nodes, n_dof)
    
    #Displacements in global coordinate system + elemental DOF displacements
    disp, compliance = displacements(K, F, dofs, n_dof)
    
    disp_local = [disp[idx] for idx in idx_expanded]
    
    F_global, f_elemental, f_axial, stress_axial, rxns = forces(n_elements, element_idx, A, k_global, disp_local, T, n_dof)
    
    #If equilibrium is not reached, exit solver
    if eq_check(rxns, F, tol) == true
        return
    end
    
    max_stress_index = findmax(abs.(stress_axial))[2]

    σ_max = stress_axial[max_stress_index]
    disp_max = findmax(disp)[1]
    
    return σ_max, disp_max, compliance
end