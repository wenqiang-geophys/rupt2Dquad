#!/bin/bash

set -e
set -x

mkdir -p data


time mpirun -np 4 ../../bin/exe_solver
