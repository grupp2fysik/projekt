#!/bin/bash
#SBATCH -A liu-compute-2026-1
#SBATCH -N 1
#SBATCH -n 4
#SBATCH -t 00:10:00
#SBATCH -J TiN_scf

module load QuantumESPRESSO/7.1

srun pw.x -in ../input/TiN_scf.in > ../output/TiN_scf.out
