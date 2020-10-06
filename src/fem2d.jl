#using Plots, LinearAlgebra
#pyplot(framestyle = :box, legend = false, color = :black)

#######################################################################################################
#######################################################################################################
#######################################################################################################
#######################################################################################################
#### Initial structures################################################################################
#######################################################################################################
#######################################################################################################
#######################################################################################################


#A node structure has an x,y position in the global space
mutable struct node
    x::Float64 #global x coordinate of node
    y::Float64 #global y coordinate of node
    dof_x::Bool #unconstrained DOF = True; fixed = False
    dof_y::Bool #unconstrained DOF = True; fixed = False
end

#An element structure has two nodes that define the start and end positions of the system
mutable struct element
    a::node #node at start
    b::node #node at end
    A::Float64 #Cross sectional area
    E::Float64 #Young's Modulus
end

#load structure indicates the node that is being loaded + global XY components of the load
mutable struct load
    n::node #node at load point
    x_val::Float64 #X component of force
    y_val::Float64 #Y component of force
end


#######################################################################################################
#######################################################################################################
#######################################################################################################
#######################################################################################################
#### Primary Functions###
#######################################################################################################
#######################################################################################################
#######################################################################################################


#Get length of element
function elem_length(elem::element)
    return sqrt((elem.b.x - elem.a.x)^2 + (elem.b.y - elem.a.y)^2)
end

#get mid point of element
function element_mid(e::element)
    x = (e.b.x - e.a.x)/2 + e.a.x
    y = (e.b.y - e.a.y)/2 + e.a.y
    return x, y
end


#Get angle of element
function θ(elem::element)
    
    #### Angle must be in degrees
    
    x_global = [1,0] #unit vector in X direciton
    vec_element = [elem.b.x - elem.a.x, elem.b.y - elem.a.y] #element as a vector
    
    #take dot product of unit x vector and the directional vector of the element
    cos = dot(x_global, vec_element) / (norm(x_global) * norm(vec_element))
    
    #get the included angle
    angle = acosd(cos)
    
    #This ensures that the clockwise-θ convention is met
    if vec_element[2] < 0
        angle = 360 - angle
    end

    return angle
end

# Which DOFs at a given node are active/constrained?
function activecheck(n::node)
    #returns tuple for DOF status of global x and y at the given node
    return (n.dof_x, n.dof_y)
end

# Single dimension array of 2n (n = number of nodes) values that indicate the active
# DOFs in order of definition in the array of nodes
function active_dofs(active_map::Array{Tuple{Bool,Bool},1})
    # active_map is an array of tuples indicating the active status of the DOF for each node
    
    #What is returned is a single 1D array of all active nodal DOFs
    #ie active_dof = [(true, true), (true, false), (false, true)]
    #gives the output [true, true, true, false, false, true]
    
    return vcat(collect.(active_map)...)
end

function node_element_link(n::node, e::element)
    #Returns a tuple of whether or not the start or end point of a given element is a given node
    return(e.a == n, e.b == n)
end


#Determine which nodes are at the start/end of a given element
function startnode(node_links::Array{Tuple{Bool,Bool},1})
    ans = findfirst([node == (true, false) for node in node_links])
    return ans
end

function endnode(node_links::Array{Tuple{Bool,Bool},1})
    ans = findfirst([node == (false, true) for node in node_links])
    return ans
end

function start_end_node_idx(element, nodes)
    nodelink = [node_element_link(node, element) for node in nodes]
    n_start = startnode(nodelink)
    n_end = endnode(nodelink)
    return [n_start, n_end]
end

# Create the vector of loads given the array of load structures
function loadvector(loads, nodes, n_dof)
    f = zeros(n_dof)
    for l in loads
        idx = 2*findfirst([node == l.n for node in nodes]) - 1
        f[idx:idx+1] = [l.x_val, l.y_val]
    end
    return f
end

## Local coordinate elemental stiffness
function k_localxy(elem::element,length::Float64)
    k = elem.A * elem.E / length .* [1 -1; -1 1]
    return k
end

## Transformation matrix
function Γ(ang)
    return [cosd(ang) sind(ang) 0 0;0 0 cosd(ang) sind(ang)]
end

function k_globalxy(k_local::Array{Float64, 2}, Γ; unit = :degrees)
    #Global_xy stiffness matrix
    k = transpose(Γ) * k_local * Γ
    
    #Note the indices of the matrix are in the form of:
    # 1: node A, global X
    # 2: node A, global Y
    # 3: node B, global X
    # 4: node B, global Y

    return k
end

## Global coordinate elemental stiffness
## This is the main function that simply outputs the global stiffness matrix
function build_K(K_dim, dof_active, nodes, elements, k_element_global)
    
    K = zeros(K_dim, K_dim)

    for i = 1:length(elements)
        
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
        k_extracted = k_element_global[i][idx_local, idx_local]
        
        K[idx_global, idx_global] .+= k_extracted
    end    
    
    #return stored matrices
    return K
end

## displacement of nodes
function d_glob(disp_long, element, nodes)
    a, b = start_end_node_idx(element, nodes)
    
    a = a*2 - 1
    b = b*2 - 1
    out_idx = [a, a+1, b, b+1]

    return disp_long[out_idx]
end

## reaction forces
function reactions(F_e, dof_list, nodes, elements)
#create empty storage matrix
F_all = zeros(length(dof_list))

#for each element
for i = 1:length(elements)
    #extract the global XY elemental end forecs
    F = F_e[i]
    
    #extract the given element information
    e = elements[i]
    
    #Determine the start and end nodes for the given element
    n_start, n_end = start_end_node_idx(e, nodes)

    #Convert to the indexes in the global DOF list 
    idx_start = 2 * n_start - 1 .+ [0,1]
    idx_end = 2 * n_end - 1 .+ [0,1]
    
    #Add the end forces of each element to the appropriate DOF:
    F_all[idx_start] .+= F[1:2]
    F_all[idx_end] .+= F[3:4]
end

#Return the force values of each DOF that are at REACTIONS
return F_all .* (dof_list .== 0)
end

## Nodes
function displaced_nodes(nodes, disp_long; scale_x = 1, scale_y = 1)
    n_disp = []
    for i = 1:length(nodes)
        n = nodes[i]
        n_new = node((n.x + disp_long[i*2 - 1]) * scale_x, (n.y + disp_long[i*2]) * scale_y, n.dof_x, n.dof_y)
        push!(n_disp, n_new)
    end
    return n_disp
end

## Elements
function displaced_elements(elements, nodes, displaced_nodes)
    elem_disp = []
    for i = 1:length(elements)
        elem = elements[i]
        a, b = start_end_node_idx(elem, nodes)
        elem_new = element(displaced_nodes[a], displaced_nodes[b], elem.A, elem.E)
        push!(elem_disp, elem_new)
    end
    return elem_disp
end

function nodecollector(nodestore::Array)
    L = length(nodestore)
    nodes = [([nodestore[i].x], [nodestore[i].y]) for i = 1:L]
    return nodes
end

function elementcollector(elementstore::Array)
    L = length(elementstore)
    elements = [([elementstore[i].a.x, elementstore[i].b.x], [elementstore[i].a.y, elementstore[i].b.y]) for i = 1:L]
    return elements
end

function loadplotter(loadstore, scale)
    if typeof(loadstore) !== Array{load,1}
        println("No loads specified/incorrect type.")
        return [([0,0], [0,0])]
    end
    #Normalize length of load arrows
    normalizer = maximum([sqrt(load.x_val^2 + load.y_val^2) for load in loadstore])
    
    L = length(loadstore)
    
    x1(i) = loadstore[i].n.x
    x2(i) = loadstore[i].n.x + loadstore[i].x_val / normalizer * scale

    y1(i) = loadstore[i].n.y
    y2(i) = loadstore[i].n.y + loadstore[i].y_val / normalizer * scale
    
    
    loads = [([x1(i), x2(i)], [y1(i), y2(i)]) for i = 1:L]
    
    return loads
end

function bc_plotter(nodes)
    #X-direction BCs
    x_fixed_idx = findall([node.dof_x == 0 for node in nodes])
    #Y-direction BCs
    y_fixed_idx = findall([node.dof_y == 0 for node in nodes])
    
    x_bcs = [([nodes[i].x], [nodes[i].y]) for i in x_fixed_idx]
    y_bcs = [([nodes[i].x], [nodes[i].y]) for i in y_fixed_idx]
    
    return x_bcs, y_bcs
end

function trussplotter(nodes, elements, loads; 
        nodelabels = false,
        node_marker = :circle,
        node_label_color = :red,
        node_label_pos = :right,
        nodesize = 3, 
        nodecolor = :black, 
        elementlabels = false,
        elementcolor = :gray,
        elementstyle = :solid,
        lwscale = 2,
        loadscale = :auto,
        bcscale = 10,
        bccolor = :green)
    
    #initiate plot
    canvas = plot()
    
    #plot Boundary Conditions
    xbc, ybc = bc_plotter(nodes)
    xmarker = [(-1, sqrt(2)/2), (0,0), (-1,-sqrt(2)/2)]
    ymarker = [(-sqrt(2)/2,-1), (0,0), (sqrt(2)/2,-1)]
    
    scatter!(xbc, marker = (Shape(xmarker), bcscale, bccolor))
    scatter!(ybc, marker = (Shape(ymarker), bcscale, bccolor))
    
    #plot elements
    #Linewidth scale
    areas = [element.A for element in elements]
    widths = transpose(areas ./ maximum(areas)) * lwscale

    plot!(elementcollector(elements), 
        color = elementcolor, 
        linestyle = elementstyle, 
        lw = widths,
        label = "")

    #plot element labels
    if nodelabels
        for i = 1:length(nodes)
            annotate!((nodes[i].x, nodes[i].y, Plots.text(string(i), node_label_color, node_label_pos)))
        end
    end

    #plot element labels
    if elementlabels
        for i = 1:length(elements)
            annotate!((element_mid(elements[i])[1], element_mid(elements[i])[2], Plots.text(string(i), :black)))
        end
    end

    #plot loads
    #Determine scale for loads
    if loadscale == :auto
        xrange = maximum(node.x for node in nodes) - minimum(node.x for node in nodes)
        yrange = maximum(node.y for node in nodes) - minimum(node.y for node in nodes)

        scale = 0.10 * max(xrange, yrange)
    else
        scale = loadscale
    end
    
    plot!(loadplotter(loads, scale), 
        color = :red, 
        arrows = true,
        label = "")
    
    #plot nodes
    scatter!(nodecollector(nodes), 
        markershape = node_marker,
        markercolor = nodecolor, 
        label = "")
    
    canvas
    return canvas
end


#######################################################################################################
#######################################################################################################
#######################################################################################################
#######################################################################################################
#### Second order Functions###
#######################################################################################################
#######################################################################################################
#######################################################################################################

function fem_init(nodes, elements)
    #extract lengths of each element
    lengths = [elem_length(e) for e in elements]
    #extract CC angle (degrees) of each element from global axis
    angles = [θ(e) for e in elements]
    
    #Determine which DOFs are active for each node (X, Y)
    dof_active = [activecheck(node) for node in nodes]
    
    #Determine array of 2*n_nodes indicating which global DOFs are active
    dof_list = active_dofs(dof_active)
    
    #Number of total DOFs (active and inactive)
    n_dof = length(dof_list)
    
    return lengths, angles, dof_active, dof_list, n_dof
end

function stiffness_init(nodes, elements, lengths, angles, dof_active, n_dof)
    #Transformation Matrix
    T = [Γ(angle) for angle in angles]
    
    #Stiffness Matrices
    k_element_local = [k_localxy(elements[i], lengths[i]) for i = 1:length(elements)]
    k_element_global = [k_globalxy(k_element_local[i], T[i]) for i = 1:length(elements)]
    
    #Global stiffness matrix
    K = build_K(n_dof, dof_active, nodes, elements, k_element_global)
    
    return T, k_element_local, k_element_global, K
end

function displacements(K, F, dof_list, n_dof)
    #Reduce K matrix and force vector to non-zero items:
    K_reduced = K[dof_list, dof_list]
    F_reduced = F[dof_list]
    
    #compute displacement vector
    disp = K_reduced \ F_reduced #K^-1 F

    compliance = transpose(F_reduced) * disp #F^T d 
    
    #Extended disp vector
    disp_long = zeros(n_dof)
    disp_long[dof_list] = disp
    
    #Convert global displacements to elemental global XY displacements
    disp_global = [d_glob(disp_long, elements[i], nodes) for i = 1:length(elements)]
    
    return disp_long, disp_global, compliance
end

function forces(nodes, elements, k_global, disp_global, T, dof_list)
    #Global element forces
    F_e = [k_global[i] * disp_global[i] for i = 1:length(elements)]
    
    #Convert global element forces t olocal element forces
    f_e = [T[i] * F_e[i] for i = 1:length(elements)]
    
    #Single axial force value
    f = [force[2] for force in f_e]
    
    #Stress values
    σ = f ./ [element.A for element in elements]
    
    #Determing Reaction Forces
    rxns = reactions(F_e, dof_list, nodes, elements)
    
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

function new_positions(nodes, elements, displacements)
    new_nodes = displaced_nodes(nodes, displacements)
    new_elements = displaced_elements(elements, nodes, new_nodes)
    
    return new_nodes, new_elements
end

function fem2d_solver(nodes, elements, loads; tol = 1e-6)
    
    #Initial geometric parameters + activity check of Degrees of Freeodm
    lengths, angles, dof_active, dof_list, n_dof = fem_init(nodes, elements)
    
    #Return global stiffness matrix (including 0 rows/columns)
    T, k_local, k_global, K = stiffness_init(nodes, elements, lengths, angles, dof_active, n_dof)
    
    #Create force vector
    F = loadvector(loads, nodes, n_dof)
    
    #Displacements in global coordinate system + elemental DOF displacements
    disp, disp_elemental, compliance = displacements(K, F, dof_list, n_dof)
    
    F_global, f_elemental, f_axial, stress_axial, rxns = forces(nodes, elements, k_global, disp_elemental, T, dof_list)
    
    #If equilibrium is not reached, exit solver
    if eq_check(rxns, F, tol) == true
        return
    end
    
    #Deformed nodes and elements
    new_nodes, new_elements = new_positions(nodes, elements, disp)
    
    return new_nodes, new_elements, disp, f_axial, stress_axial, rxns, compliance
end

function fem2d_solver_modified(nodes, elements, loads; tol = 1e-6)
    
    #Initial geometric parameters + activity check of Degrees of Freeodm
    lengths, angles, dof_active, dof_list, n_dof = fem_init(nodes, elements)
    
    #Return global stiffness matrix (including 0 rows/columns)
    T, k_local, k_global, K = stiffness_init(nodes, elements, lengths, angles, dof_active, n_dof)
    
    #Create force vector
    F = loadvector(loads, nodes, n_dof)
    
    #Displacements in global coordinate system + elemental DOF displacements
    disp, disp_elemental, compliance = displacements(K, F, dof_list, n_dof)
    
    return compliance
end



##############################################################################
##############################################################################
#################GRADIENT METHODS#############################################

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
        

# ##############################################################################
# ##############################################################################
# #################COMPLIANCE OPTIMIZER#########################################

# function opt_compliance(nodes, elements, loads, area_init, V_max; 
#     write_sol = false,
#     x_lowerbound = 0.0,
#     xtol = 1e-6, 
#     ftol = 1e-4)

#     n_elements = length(elements)

#     #Initialize optimization model
#     compliance_opt = Opt(:LD_MMA, n_elements)
#     compliance_opt.lower_bounds = x_lowerbound * ones(n_elements)
#     compliance_opt.xtol_rel = xtol
#     compliance_opt.ftol_abs = ftol

#     #sub function for the objective function
#     function compliance_areas(ars, grad)
#         new_elems = [element(elements[i].a, elements[i].b, ars[i], elements[i].E) for i = 1:length(ars)]
#         compliance = fem2d_solver_modified(nodes, new_elems, loads)
#         grad_store = adjoint_gradient(nodes, new_elems, loads)
        
#         if length(grad) > 0
#             for i = 1:length(grad)
#                 grad[i] = grad_store[i]
#             end
#         end
#         return compliance
#     end

#     #add objective function to model
#     compliance_opt.min_objective = compliance_areas

#     #Volume constraint
#     function vol_constraint(ars, grad)
#         new_elems = [element(elements[i].a, elements[i].b, ars[i], elements[i].E) for i = 1:length(ars)]
#         grad_store = adjoint_gradient(nodes, new_elems, loads)
        
#         if length(grad) > 0
#             for i = 1:length(grad)
#                 grad[i] = grad_store[i]
#             end
#         end
#         return sum(lengths .* ars) * 12  - V_max
#     end

#     #Add constraint to model
#     inequality_constraint!(compliance_opt, vol_constraint, 1e-6)

#     #solve
#     (minf,minx,ret) = optimize(compliance_opt, area_init)

#     if write_sol
#         numevals = compliance_opt.numevals # the number of function evaluations
#         println("got $minf at $minx after $numevals iterations (returned $ret)")
#     end

#     return minf, minx, ret
# end
