#!/bin/bash

set -e
set -x

fnm=bend
gmsh -2 ${fnm}.geo -o ${fnm}.inp
meshio convert ${fnm}.inp ${fnm}.exo
rm -f ${fnm}.inp
