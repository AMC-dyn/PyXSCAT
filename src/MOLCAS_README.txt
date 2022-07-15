Hello!

There are two things required to use a molcas calculations - the output file (e.g. molcas.log) and a molden file (e.g. molcas.rasscf.molden)

From the former we read the ci vector, and the from the latter we read the orbitals/basis set/geometry.

These both have special readers - molcas_ci_reader.py and molden_reader_nikola_molcas.py.

The interface can currently read both the CASSCF/RASSCF output (of the program `rasscf`) and the CASPT2 output (of the eponymous program). The CASPT2 orbitals are *without* the PMCAS perturbative correction, and so represent an approximation of the total wavefunction in the reference space, i.e. the set of states generated in the CAS/RAS calculation.

The program should automatically detect the symmetry, but does not run using the symmetry - it takes the symmetric wavefunctions/ci vectors and runs as if they were not symmetric. This still represents a computational speed up as the total space of configurations/determinants is reduced.

**IMPORTANT**
You must use `PRWF = -1.` in your input scripts, or you will limit the total space of determinants. The code currently cannot handle an incorrect number of determinants, but this could be changed fairly easily.

The interface converts from the CSF output given by molcas to a SD basis. This might take some time (especially if lots of CSFs are used).


