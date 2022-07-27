#!/bin/bash -l

#SBATCH -n 1
#SBATCH -t 3-00:00:00
#SBATCH --cpus-per-task=8

conda activate pymodsim

pymodsim -n 23 -s scnba_human.fasta -p AF-Q9UI33-F1-model_v2.pdb -t in
