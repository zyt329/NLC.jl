# NLC

[![Build Status](https://github.com/zyt329/NLC.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/zyt329/NLC.jl/actions/workflows/CI.yml?query=branch%3Amain)

This package intends to provide with a collection of functions that allows one do the Numerical Lined Cluster (NLC) expansion painlessly. 

# Version

This version highlights the following features:

[1] This version works with the '.json' cluster information provided by Pranav

[2] This version stores diagonalization data for later use. The stored data were only E and M. With compression, the data should take up ~20GB to the 14th order for triangular lattice Heisenberg model.

[3] To check: This version uses not only Sz symmetry, but also parity symmetry, for the Heisenberg model. This reduces the computation by a little more than a factor of 2. 

[4] Implemented efficient parallel computing with MPI for different clusters. Now there should be minimal conflict in diagonalizing large number of clusters.

# Installation Instructions

Here's instruction on how to install the package and run the code:

[1] git clone the package to wherever you like on your machine.

[2] install julia on your machine

[3] navigate your terminal to the package's folder (PathToYourFolder/NLC.jl). 

[4] open julia repl. Then press the "]" button, entering the package manager.

[5] type "dev ." (without the quotation marks) and hit Enter.

[6] Once the last command finishes running, press the Backspace button, going back to the normal julia repl. Job Done! The package should have been installed properly now.

# To run some examples:

Navigate into the 'examples' folder. The example files allow you to run NLC simulation of the XXZ model on a triangular lattice to the 9th order. The procedure to run a full NLC simulation consists of 4 steps, and is accomplished by running 4 files:

[1] Diagonalize the clusters: run 'julia xxz_diagonalization.jl'

This would automatically create a folder and store the diagonalization data within. Some example cluster information are stored in the 'cluster_info' folder.

[2] Thermalize the clusters with diagonalization data: run 'julia xxz_simulation.jl'

This takes thermal averages of all clusters for the quantities to be measured.

[3] Do the NLC sum with data in [2]: run 'julia xxz_NLC_sum.jl'

[4] Do resummation to improve convergence: run 'julia xxz_resum.jl'

The script would do the Wynn-epsilon and Euler resummation to all orders possible. Raw sums (result without resummation) will also be stored after running resummation.

The final result after resummation can be accessed in the folder 'test_simulation_J_z[J_xy1.0000-1/resummation-1/'


# Package Dependencies

All packages except the following 2 should be installed once NLC.jl is installed on the machine following the installation instructions. The following 2 need to be installed separately:

[1] JLD2

This package is used to store and read data from disk in this package.

[2] CodecZlib

This also needs to be installed separately. The package is used by JLD2.jl to compress diagonalization data. Somehow it's not in the dependency of JLD2 and needs to be installed separately.