import pandas as pd
import numpy as np
import molden_reader_nikola_morder as mldr
from scipy.io import FortranFile
from copy import deepcopy
i = 0
nmomax =7

gtos, atoms, coeffs, mos, groupC,ncontr = mldr.read_orbitals('lif_ac_m.mld', N=nmomax, decontract=False)


print('MO shape: ', mos.shape)

mmodt = np.transpose(mos)

trdm = np.zeros((nmomax, nmomax, nmomax, nmomax))
with open(file='twordm_fortran_bit_2.dat') as f:
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



print('trdm',trdm[0,0,0,0])
mmod=mos

print(mmod[0,:])
CC = np.einsum("up,pqrs->uqrs", mmod, trdm, optimize='optimal')
print('number to compare 1 : ',CC[1,2,2,6])


CC.tofile('prueba2.dat')
# f = FortranFile('prueba2.dat', 'w')
# f.write_record(CC)
# f.close()
print(CC.shape)
CC = np.einsum("vq,uqrs->uvrs", mmod, CC, optimize='optimal')
print('number to compare 1 : ',CC[1,2,2,6])
print(CC.shape)
CC = np.einsum("kr,uvrs->uvks", mmod, CC, optimize='optimal')
print('number to compare 1 : ',CC[1,2,2,6])
print(CC.shape)
CC1 = np.einsum("ls,uvks->uvkl", mmod, CC, optimize='optimal')
print('number to compare 1 : ',CC1[1,2,2,6])
#CC1.tofile('Zcotr.dat')
print('number to compare 1: ',CC1[0,1,1,2])
print('number to compare 2: ',CC1[1,1,2,0])
print('number to compare 3: ',CC1[1,0,3,2])
print('number to compare 4: ',CC1[2,3,1,0])
print('number to compare 5: ',CC1[2,3,0,1])
print('number to compare 6: ',CC1[0,1,3,2])
print('number to compare 7: ',CC1[3,2,0,1])
print('number to compare 8: ',CC1[3,2,1,0])
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
print(Zbig[1,2,2,6])
Zbig.tofile('Zcotr.dat')
print(Zbig[0,1,2,3])


