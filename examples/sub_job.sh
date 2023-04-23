#!/bin/bash
#SBATCH --job-name=NLC_xxz 
#SBATCH --output=NLC_xxz.log
#SBATCH --partition=puma-i9
for i in 1 2 4 8 16
do
   time julia xxz_simulation.jl $i > runs.out 
done

