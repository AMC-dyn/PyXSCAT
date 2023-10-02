import pandas as pd
import numpy as np
import molden_reader_nikola as mldr
from scipy.io import FortranFile
from copy import deepcopy
i = 0
nmomax = 18

gtos, atoms, coeffs, mos, groupC = mldr.read_orbitals('co_cas.mld', N=nmomax, decontract=False)


print('MO shape: ', mos.shape)

mmodt = np.transpose(mos)

trdm = np.zeros((nmomax, nmomax, nmomax, nmomax))
with open(file='twordm_fortran.dat') as f:
    lines = f.readlines()
    print(len(lines))
    totalfin = np.zeros(len(lines))
    nlines = len(lines)
    m1 = np.ones(nlines)
    m2 = np.ones(nlines)
    m3 = np.ones(nlines)
    m4 = np.ones(nlines)

    for line in lines:
        all_num = line.strip().split()
        m1[i] = int(all_num[0]) - int(1)
        m2[i] = int(all_num[1]) - int(1)
        m3[i] = int(all_num[2]) - int(1)
        m4[i] = int(all_num[3]) - int(1)
        totalfin[i] = np.float64(all_num[4])
        trdm[int(m1[i]), int(m2[i]), int(m3[i]), int(m4[i])] = totalfin[i]
        i += 1
print('trying a masiva calculacion')




mmod=mos


CC = np.einsum("up,pqrs->uqrs", mmod, trdm, optimize='optimal')
print(CC.shape)
CC = np.einsum("vq,uqrs->uvrs", mmod, CC, optimize='optimal')
print(CC.shape)
CC = np.einsum("kr,uvrs->uvks", mmod, CC, optimize='optimal')
print(CC.shape)
CC1 = np.einsum("ls,uvks->uvkl", mmod, CC, optimize='optimal')
#CC1.tofile('Zcotr.dat')

CC = CC1

CC = CC + np.einsum('abcd->bacd', CC1, optimize='optimal')

CC = CC + np.einsum('abcd->abdc', CC1, optimize='optimal')

print('up to here')

CC = CC + np.einsum('abcd->badc', CC1, optimize='optimal')

CC = CC + np.einsum('abcd->cdab', CC1, optimize='optimal')

print('half way')
CC = CC + np.einsum('abcd->dcab', CC1, optimize='optimal')

CC = CC + np.einsum('abcd->cdba', CC1, optimize='optimal')

CC = CC + np.einsum('abcd->dcba', CC1, optimize='optimal')

Zbig = CC / 8.0000
Zbig.tofile('Zcotr.dat')
print(CC1.shape)

