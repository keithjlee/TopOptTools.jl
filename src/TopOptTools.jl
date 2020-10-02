module TopOptTools
using Plots, LinearAlgebra, ForwardDiff, DataFrames
pyplot(framestyle = :box, legend = false, color = :black)

include("fem2d.jl")
include("opt.jl")

#Export data types
export node, element, load

#Export main solver functions
export fem2d_solver, fem2d_solver_modified

#Export Topology optimization method
export minimize_unconstrained

end
