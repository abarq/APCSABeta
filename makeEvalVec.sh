#!/bin/bash

FLAGS=("-fno-tree-vectorize" "-ftree-vectorize -fopt-info-vec-optimized" "-ftree-vectorize -fopt-info-vec-optimized" "-ftree-vectorize -march=native -mtune=native -fopt-info-vec-optimized")
casesCostFunction=(3 2 1 1)
exeName="evalVectorization"

n=1

for caso in "${FLAGS[@]}"
do
gfortran -O3  $caso  -o $exeName kinds.f90 randomNum.f90 inputData.f90 instanceProblem.f90 APCSA.f90 evalVectorization.f90
./$exeName ${casesCostFunction[$(($n-1))]}  $n
#echo ${casesCostFunction[$(($n-1))]}
n=$(($n+1))


done
