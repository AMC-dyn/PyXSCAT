import pandas as pd
import numpy as np
import molden_reader_nikola as mldr
from scipy.io import FortranFile
from copy import deepcopy
i = 0
nmomax = 10

gtos, atoms, coeffs, mos, groupC = mldr.read_orbitals('co_cas.mld', N=nmomax, decontract=False)
mmod = mos
print('MO shape: ', mos.shape)

mmodt = np.transpose(mos)
print(groupC)

print(coeffs)
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

i = 0

maxorb = np.max(m1)
# rata4.tofile('Zcotr.dat')
m1 = np.asarray(m1, dtype=np.int32)
m2 = np.asarray(m2, dtype=np.int32)
m3 = np.asarray(m3, dtype=np.int32)
m4 = np.asarray(m4, dtype=np.int32)
ii = 0
jj = 1
k = 1
r = 0
mmodt = np.transpose(mmodt)

mmod = np.transpose(mmod)

totalfin = totalfin
print('prueba:', np.sum(totalfin * ((mmod[m1, ii] * mmod[m2, jj]) * (mmod[m3, k] * mmod[m4, r]))))
temp1 = totalfin * ((mmod[m1, ii] * mmod[m2, jj]) + (mmod[m1, jj] * mmod[m2, ii]))
temp2 = (mmod[m3, k] * mmod[m4, r]) + (mmod[m3, r] * mmod[m4, k])

temp3 = totalfin * (mmod[m1, r] * mmod[m2, k] + mmod[m1, k] * mmod[m2, r])
temp4 = (mmod[m3, ii] * mmod[m4, jj]) + (mmod[m3, jj] * mmod[m4, ii])
mmod = np.transpose(mmod)
Zcont1 = np.sum(temp3 * temp4 + temp1 * temp2)
print('resultado: ', Zcont1)

CC = np.einsum("up,pqrs->uqrs", mmod, trdm, optimize='greedy')
print(CC.shape)
CC = np.einsum("vq,uqrs->uvrs", mmod, CC, optimize='greedy')
print(CC.shape)
CC = np.einsum("kr,uvrs->uvks", mmod, CC, optimize='greedy')
print(CC.shape)
CC1 = np.einsum("ls,uvks->uvkl", mmod, CC, optimize='greedy')
# CC1.tofile('Zcotr.dat')

# CC1.tofile('Zcotr.dat')
del CC, temp1, temp2, temp3, temp4, m1, m2, m3, m4, totalfin, mmod, mmodt
print(CC1.shape)
print('intermediate step')
print('prueba1 :', CC1[1, 1, 1, 1] * coeffs[0] * coeffs[1] * coeffs[1] * coeffs[0])
# print((CC1[7, 6, 6, 10] + CC1[7, 6, 10, 6] + CC1[6, 7, 6, 10] + CC1[6, 7, 20, 6] + CC1[20, 6, 6, 7] + CC1[6, 20, 6, 7] + CC1[
#  20, 6, 7, 6] + CC1[6, 20, 7, 6]) / 8.0)


print('CC1:', CC1[4, 0, 0, 0] * coeffs[0] * coeffs[1] * coeffs[1] * coeffs[16:20])
print(CC1[0, 1, 1, 2])
print(CC1[1, 0, 1, 2])
print(CC1[0, 1, 2, 1])

CC = CC1

CC = CC + np.einsum('abcd->bacd', CC1, optimize='greedy')

CC = CC + np.einsum('abcd->abdc', CC1, optimize='greedy')

print('up to here')

CC = CC + np.einsum('abcd->badc', CC1, optimize='greedy')

CC = CC + np.einsum('abcd->cdab', CC1, optimize='greedy')

print('half way')
CC = CC + np.einsum('abcd->dcab', CC1, optimize='greedy')

CC = CC + np.einsum('abcd->cdba', CC1, optimize='greedy')

CC = CC + np.einsum('abcd->dcba', CC1, optimize='greedy')

Zbig = CC / 8.0000
print(CC.shape)
rr =deepcopy(coeffs)
# print(coeffs)
rr[10] = coeffs[10]
rr[11] = coeffs[13]
rr[12] = coeffs[16]
rr[13] = coeffs[11]
rr[14] = coeffs[14]
rr[15] = coeffs[17]
rr[16] = coeffs[12]
rr[17] = coeffs[15]
coeffs = deepcopy(rr)
print(coeffs)
print('Zbig:', Zbig[0, 0, 0, 5] * coeffs[0] * coeffs[1] * coeffs[1] * coeffs[16:19])
Zbig.tofile('Zcotr.dat')


print(CC[0, 1, 1, 0])
print(CC[1, 2, 2, 1])

mmod = gtos.mo
print(mmod[:,0])
print(mmod[12,0])
rr =deepcopy(mmod)
print(rr[:,0]-mmod[:,0])

rr[10, :] = mmod[10, :]
rr[11, :] = mmod[13, :]
rr[12, :] = mmod[16, :]
rr[13, :] = mmod[11, :]
rr[14, :] = mmod[14, :]
rr[15, :] = mmod[17, :]
print(rr[16,0],mmod[12,0])
rr[16, :] = mmod[12, :]
print(rr[16,0],mmod[12,0])
rr[17, :] = mmod[15, :]
rr[18, :] = mmod[18, :]

mmod =deepcopy(rr)

print(mmod[:,0])

print(mmod.shape)
CC = np.einsum("up,pqrs->uqrs", mmod, trdm, optimize='greedy')
print(CC.shape)
CC = np.einsum("vq,uqrs->uvrs", mmod, CC, optimize='greedy')
print(CC.shape)
CC = np.einsum("kr,uvrs->uvks", mmod, CC, optimize='greedy')
print(CC.shape)
CC1 = np.einsum("ls,uvks->uvkl", mmod, CC, optimize='greedy')

CC = CC1

CC = CC + np.einsum('abcd->bacd', CC1, optimize='greedy')

CC = CC + np.einsum('abcd->abdc', CC1, optimize='greedy')

print('up to here')

CC = CC + np.einsum('abcd->badc', CC1, optimize='greedy')

CC = CC + np.einsum('abcd->cdab', CC1, optimize='greedy')

print('half way')
CC = CC + np.einsum('abcd->dcab', CC1, optimize='greedy')

CC = CC + np.einsum('abcd->cdba', CC1, optimize='greedy')

CC = CC + np.einsum('abcd->dcba', CC1, optimize='greedy')

Zbig = CC / 8.0000
print(CC1.shape)
print(Zbig[0, 1, 1, 10:20])
