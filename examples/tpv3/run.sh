#!/bin/bash

set -e
set -x

rm -rf data/fault*
rm -rf data/wave*

time mpirun -np 4 ../../bin/exe_solver
