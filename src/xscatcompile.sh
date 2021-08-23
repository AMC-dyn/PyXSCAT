#!/bin/bash

rm *.o *.mod
gfortran -lblas -O3 -c maricones.f95
f2py3 -lblas -c --fcompiler=gfortran -I. maricones.o Integrals_ijkr.f95 -m integrals_wrapper only: total_scattering_calculation
