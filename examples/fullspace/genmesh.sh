#!/bin/bash

set -e
set -x

fnm=rect
gmsh -2 ${fnm}.geo -o ${fnm}.inp
meshio convert ${fnm}.inp ${fnm}.exo
rm -f ${fnm}.inp
