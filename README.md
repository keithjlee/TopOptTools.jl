# TopOptTools
This is a series of tools developed from scratch in Julia both as a learning method as well as functional use for the MIT class 1.583: Topology Optimization of Structures taken in Fall 2020. The main goal of this repository is for me to learn best-practices for coding, both in the language itself (taking full advantage of multiple dispatch and data structures), as well as properly documenting, version-controlling, and maintaining a code base.

The base package of `TopOptTools.jl` contains increasing sub-packages for finite element analysis and topology optimization tools, mostly for 2D trusses. A lot of this is hacked together, and I know there are much more elegant (and efficient) ways to approach node-element connectivity structures, so much of this will hopefully change.

The intent is to slowly develop these tools into more robust, general systems that can be immediately reused for my own research in structural design.

## A quick guide on using `opt_unconstrained()`
`opt_unconstrained` is the unconstrained gradient-based optimizer. The required inputs are:

1. `func`: This is the input of type *function*. 
   1. Any Julia-readable function can be used.
2. `x_start`: This is an array whose length equals the number of required inputs in `func`.
   1. There is currently no failsafe if these arrays do not line up, and the error will occur elsewhere
3. `descent_method`: This is a symbolic input with three options:
   1. `:steepest`: Steepest descent method
   2. `:cg`: Conjugate gradient method
   3. `:newton`: Newton's method
4. `iter_max`: this is an Int value that indicates the maximum number of loops allowed
5. `alpha_max`: this is the maximum step size allowed in the *line-search* algorithm, where currently only the reducing interval method is implemented
6. `epsilon`: this is the allowable function error for optimality

One optional input, `annotate`, which is default set to `true`, prints out step-by-step results of the solver. 

## Small example:
Say we want to find the unconstrained minimum of: $(\vec{x}) = 25x_1^3 + sqrt(x_2) - cos(x_3pi)$,

Starting at the point $\vec{x} = [1.2, 3.0, 6.4]^T$.

1. First define the express the objective function as a Julia function:
   1. `f(x) = 25x[1]^3 + sqrt(x[2]) - cos(x[3] * pi)`
2. Then define the starting point as a vector:
   1. `x_0 = [1.2, 3.0, 6.4]`
3. Then input into the optimization function with the select descent and tolerance limits:
   1. `solution = opt_unconstrained(f, x_0, :cg, 10_000, 1, 0.001; annotate = false)`

The function will return a `DataFrame` type with the following information:
1. Total number of iterations
2. The value of the design variables at the minimizer
3. The value of the gradient at each iteration
4. The value of the norm of the gradient at each iteration
5. The directional vector at each iteration
6. The step size at each iteration
7. The number of iterations in the line-search subroutine