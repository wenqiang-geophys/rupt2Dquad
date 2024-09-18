#!/bin/bash

set -e
set -x

rm -rf output && mkdir output

time mpirun -np 4 ../../bin/exe_solver
