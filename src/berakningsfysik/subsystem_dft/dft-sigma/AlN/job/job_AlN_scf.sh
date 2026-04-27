#!/bin/bash
#SBATCH -A liu-compute-2026-1
#SBATCH -N 1
#SBATCH -n 4
#SBATCH -t 00:10:00
#SBATCH -J AlN_scf

module load QuantumESPRESSO/7.1

srun pw.x -in ../input/AlN_scf.in > ../output/AlN_scf.out
