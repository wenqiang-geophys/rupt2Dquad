#!/bin/bash

set -e
set -x

fnm=tpv5_2d_symm
fnm=tpv5_2d
gmsh -2 ${fnm}.geo -o ${fnm}.inp
meshio convert ${fnm}.inp ${fnm}.exo
#rm -f ${fnm}.inp
