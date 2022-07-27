#!/bin/bash -l

#SBATCH -n 1
#SBATCH -t 3-00:00:00
#SBATCH --cpus-per-task=8

conda activate pymodsim

pymodsim -n 23 -s scnba_human.fasta -p AF-Q9UI33-F1-model_v2.pdb -N 65 -C 1613 -l 418-554,811-1032,1319-1341 -f 3.0 -t in
