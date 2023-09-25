MAKEFLAGS += --no-builtin-rules --no-buliltin-variables

name = Scattering_code

FC = ifort 

FCFLAGS = -O3 -xhost -fPIC -align array64byte  -qopenmp -heap-arrays  -qmkl

COMPILE.f08 = $(FC) $(FCFLAGS) -c

MAINFILE=Main.f90

SOURCES= calculate_form_factors.f90 unique_module.f90 twordms.f90 one_rdm.f90 Integrals.f90 Variables.f90 main_calculation.f90 linspace.f90

mainprog: $(subst .f90,.o,$(SOURCES))
	$(FC) $(FCFLAGS) $(MAINFILE) -o $@ $+

.PHONY: clean
clean:
	-rm -f *.o *.mod *.smod main

%.o:%.f90 
	$(COMPILE.f08) -o $*.o $< 
	@touch $@


