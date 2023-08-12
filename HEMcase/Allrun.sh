#!/bin/bash

source /media/gagan/OpenFOAM/openfoam-OpenFOAM-v2006/etc/bashrc

foamListTimes -rm

rm -rf 0

blockMesh

cp -r 0.orig 0

decomposePar -force

foamJob -parallel -screen  rhoCentralFoamHEM
