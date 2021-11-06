#!/bin/bash

#SBATCH --partition=full
#SBATCH --job-name=IMT2112

#SBATCH --output=log.out
#SBATCH --ntasks=3
#SBATCH --cpus-per-task=1

mpic++ T3.cpp
mpirun b.out
