#!/bin/bash
#SBATCH --job-name=tpv14_test
#SBATCH --output=slurm_gp229_%j.out
#SBATCH --ntasks=4
#SBATCH --ntasks-per-node=4
#SBATCH --cpus-per-task=1
#SBATCH --time=00:10:00
#SBATCH --partition=serc

ml flexiblas netcdf-fortran

set -e
set -x

rm -rf output && mkdir output

srun ../../bin/exe_solver

