#!/bin/sh
gfortran -c calculate_form_factors.f90
gfortran -c besseldiff1.f95
gfortran test.f95 calculate_form_factors.o besseldiff1.o -o diff1test.exe
