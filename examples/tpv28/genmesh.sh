#!/bin/bash

set -e
set -x

fnm=tpv28_2d
gmsh -2 ${fnm}.geo -o ${fnm}.inp
meshio convert ${fnm}.inp ${fnm}.exo
rm -f ${fnm}.inp
