#!/bin/bash
module load netcdf-fortran
module load flexiblas
mkdir -p bin obj
make clean
make -f Makefile.sherlock
