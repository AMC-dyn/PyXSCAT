#!/bin/sh
gfortran -c calculate_form_factors.f90
gfortran -c besseldiff2.f95
gfortran test.f95 calculate_form_factors.o besseldiff2.o -o diff2test.exe
