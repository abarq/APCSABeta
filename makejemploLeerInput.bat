@ECHO OFF

gfortran -O3  -ftree-vectorize -march=native -mtune=native -fopt-info-vec-optimized  -o ejecutable kinds.f90 randomNum.f90 inputData.f90 instanceProblem.f90 APCSA.f90 parametersTuning.f90 

pause