#!/bin/bash
#SBATCH --job-name=test1
#SBATCH --output=slurm_gp229_%j.out
#SBATCH --ntasks=4
#SBATCH --ntasks-per-node=4
#SBATCH --mem=1GB
#SBATCH --cpus-per-task=1
#SBATCH --time=00:10:00
#SBATCH --partition=serc

ml flexiblas netcdf-fortran

set -e
set -x

rm -f data/wave*
rm -f data/fault*
rm -f data/body*

srun ../../bin/exe_solver

