#!/bin/bash

set -e
set -x


time mpirun -np 4 ../../bin/exe_solver
