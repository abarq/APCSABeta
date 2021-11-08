#!/bin/bash

FLAGS=( "-ftree-vectorize -march=native   -mtune=native -fopt-info-vec-optimized")
exeName="evalVectorization"

gfortran -O3  $FLAGS -pg  -o $exeName kinds.f90 randomNum.f90 inputData.f90 instanceProblem.f90 APCSA.f90 evalVectorization.f90




