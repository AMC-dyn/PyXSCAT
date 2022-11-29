
#/bin/bash 

f1='mpiifort -O3 -fp-model fast=2 -xCORE-AVX512 -align array64byte -heap-arrays  -qmkl  -c'

rm -rf *.mod *.o

$f1 calculate_form_factors.f90
$f1 unique_module.f90
$f1 twordms.f90
$f1 one_rdm.f90
$f1 Integrals_mpi.f90
$f1 Variables.f90

$f1 main_calculation.f90 
$f1 linspace.f90
mpiifort -O3  -fp-model fast=2 -xCORE-AVX512 -align array64byte -heap-arrays -qmkl Main_mpi.f90 *.o -o  Main.exe

