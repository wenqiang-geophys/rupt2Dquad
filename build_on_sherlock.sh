#!/bin/bash
module load netcdf-fortran
module load flexiblas
mkdir -p bin obj
make -f Makefile.sherlock
