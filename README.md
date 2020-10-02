# TopOptTools
This is a series of tools developed from scratch in Julia both as a learning method as well as functional use for the MIT class 1.583: Topology Optimization of Structures taken in Fall 2020. The main goal of this repository is for me to learn best-practices for coding, both in the language itself (taking full advantage of multiple dispatch and data structures), as well as properly documenting, version-controlling, and maintaining a code base.

The base package of `TopOptTools.jl` contains increasing sub-packages for finite element analysis and topology optimization tools, mostly for 2D trusses. A lot of this is hacked together, and I know there are much more elegant (and efficient) ways to approach node-element connectivity structures, so much of this will hopefully change.

The intent is to slowly develop these tools into more robust, general systems that can be immediately reused for my own research in structural design.

## A quick guide on using `opt`
`