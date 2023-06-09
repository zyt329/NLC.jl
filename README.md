# NLC

[![Build Status](https://github.com/zyt329/NLC.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/zyt329/NLC.jl/actions/workflows/CI.yml?query=branch%3Amain)

This package intends to provide with a collection of functions that allows one do the Numerical Lined Cluster (NLC) expansion painlessly. 

# Version

This version highlights the following features:

[1] This version works with the '.json' cluster information provided by Pranav

[2] This version stores diagonalization data for later use. The stored data were only E and M. With compression, the data should take up ~20GB to the 14th order.

*[3] To check: This version uses not only Sz symmetry, but also parity symmetry. This reduces the computation by a little more than a factor of 2. Whether the implementation have bug needs to be checked. 

*[4] To check: I have incorporate parallel computing for the diagonalization in a naive way that the program will diagonalize a cluster only if it can't find the data of the cluster in the folder. The reliability of this implementation needs to be checked.