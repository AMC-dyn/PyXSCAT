
#/bin/bash 

f1='ifort -O3  -fPIC -fopenmp -heap-arrays -mkl -c'

rm -rf *.mod *.o

$f1 unique_module.f90
$f1 twordms.f90
$f1 one_rdm.f90
$f1 Integrals.f90
$f1 Variables.f90

f2py3 -c --fcompiler=ifort --f90flags='-qopenmp -mkl -O3  -heap-arrays -std=f2003' -I.  *.o main_calculation.f90 -m integrals_wrapper

