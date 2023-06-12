# NLC

[![Build Status](https://github.com/zyt329/NLC.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/zyt329/NLC.jl/actions/workflows/CI.yml?query=branch%3Amain)

This package intends to provide with a collection of functions that allows one do the Numerical Lined Cluster (NLC) expansion painlessly. 

# Version

This version highlights two main features (comparing to other versions):

[1] This version works with the '.json' cluster information provided by Pranav

[2] This version incorporates diagonalization into the Thermal averaging process. This means the diagonalization information will not be stored, so one would have to run diagonalization each time of the simulations. However, this does save up the disk space which could be a problem for NLC calculations to the 14th order on a triangular lattice. (whose diagonalization data would by my estimation take up ~200GB)

*[3] To do: incorporate parallel computing for the diagonalization process.  