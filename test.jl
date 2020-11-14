include("optv2.jl")

nodes = [[0,0], [4,4], [8,0], [12,4]]
eidx = [[1,2], [1,3], [2,3], [2,4], [3,4]]
dofs = [false, false, true, true, true, false, true, true]

loads = [[0, -1], [0, -1]]
positions = [[4,4], [12,4]]

A = 10 * ones(length(eidx))
E = 20 * ones(length(eidx))

stress, d, comp = analysis(nodes, dofs, eidx, A, E, loads, positions)
