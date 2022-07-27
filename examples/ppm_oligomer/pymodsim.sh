#!/bin/bash -l

#SBATCH -n 1
#SBATCH -t 3-00:00:00
#SBATCH --cpus-per-task=8

conda activate pymodsim

pymodsim -n 3 -p 10mer_prot_base.pdb -t out -c A,B,C,D
