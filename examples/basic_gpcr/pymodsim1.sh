#!/bin/bash -l

#SBATCH -n 1
#SBATCH -t 3-00:00:00
#SBATCH --cpus-per-task=8

conda activate pymodsim

pymodsim -s aa2br_human.fasta -t out
