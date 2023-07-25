
#/bin/bash 

f1='ifort -O2 -align array64byte  -xCORE-AVX512 -qopenmp  -qmkl -heap-arrays -c'

rm -rf *.mod *.o

$f1 calculate_form_factors.f90
$f1 unique_module.f90
$f1 twordms.f90
$f1 one_rdm.f90
$f1 Integrals.f90
$f1 Variables.f90

$f1 main_calculation.f90 
$f1 linspace.f90
 ifort -O2 -align array64byte  -qopenmp -xCORE-AVX512  -qmkl -heap-arrays Main.f90 *.o -o  Main.exe
