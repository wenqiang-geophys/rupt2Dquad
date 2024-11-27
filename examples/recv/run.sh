#!/bin/bash

set -e
set -x

rm -f data/wave*
rm -f data/fault*

time mpirun -np 4 ../../bin/exe_solver
