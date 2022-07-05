
#/bin/bash 

f1='ifort -O3 -fast -fPIC -fopenmp -heap-arrays -unroll=3 -xHost -ipo -qmkl  -c'

rm -rf *.mod *.o
$f1 calculate_form_factors.f90
$f1 unique_module.f90
$f1 twordms.f90
$f1 one_rdm.f90
$f1 Integrals.f90
$f1 Variables.f90

f2py3 -c --fcompiler=ifort --f90flags='-qopenmp -fast -qmkl -O3 -unroll=3 -xHost -ipo -heap-arrays -std=f2003' -I.  *.o main_calculation.f90 -m integrals_wrapper

