#!/bin/sh
gfortran -c calculate_form_factors.f90
gfortran -c besseldiff0.f95
gfortran test.f95 calculate_form_factors.o besseldiff0.o -o diff1test.exe
