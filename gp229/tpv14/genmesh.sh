#!/bin/bash

set -e
set -x

fnm=tpv14_2d
gmsh -2 ${fnm}.geo -o ${fnm}.inp
meshio convert ${fnm}.inp ${fnm}.exo
rm -f ${fnm}.inp
