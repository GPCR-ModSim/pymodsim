#!/bin/bash -l

#SBATCH -n 1
#SBATCH -t 3-00:00:00
#SBATCH --cpus-per-task=8

conda activate pymodsim

pymodsim -n 23 -s aa2br_human.fasta -p ranked_0.pdb -l 0 -t out