
#/bin/bash 

f1='ifort -O3 -xhost -fPIC -align array64byte  -qopenmp -heap-arrays  -qmkl -c'

#rm -rf *.mod *.o

$f1 Reader.f90
$f1 calculate_form_factors.f90
$f1 unique_module.f90
$f1 twordms.f90
$f1 one_rdm.f90
#$f1 Integrals.f90
#$f1 Variables.f90
$f1 p0_cases.f90
$f1 MD.f90
$f1 Bessels_j0.f90
$f1 Zcontr.f90
#$f1 Total_scattering_j0_contr.f90
#$f1 Total_scattering_j0_groups.f90
$f1 Total_scattering_j0_groups_fast.f90
$f1 Total_electron_scattering_fast.f90
#$f1 main_calculation.f90 
$f1 linspace.f90
 ifort -O3 -align array64byte  -xhost -fPIC -qopenmp  -heap-arrays  -qmkl  Main.f90 *.o -o  Main2.exe
