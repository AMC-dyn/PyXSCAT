import numpy as np
import molden_reader_nikola_morder as mldreader
import molden_reader_nikola_molcas as mcmldreader
import molden_reader_nikola_bagel as bgmldreader
import molcas_ci_reader as mc
import twordm_red as td
from textwrap import wrap

# If there is an external file with CIvecs or 2rdms JeremyR==True
molpro = True
bagel = False
molcas = False
extra_precision_molcas = True
caspt2 = False
jeremyR = False
fileJeremy = 'QC-2RDM-AD.dat'
# If we convert form a MCCI calculation to a bitwise operation mcci==True
mcci = False
# If we have a HF calculation and our Civector is represented by a single determinant hf==True
hf = False
# States involved
state1 = 1
state2 = 1
closed = 2
qmin = 1E-10
qmax = 10
npoints = 100
cutoffcentre = 0.1  # suggested: cutoffcentre = 0.01;
# the cutoff for the Z integral
cutoffz = 1e-20  # suggested: cutoffz = 1E-9;
# the cutoff for the product of the MD coefficients
cutoffmd = 1e-20  # suggested: cutoffmd = 1E-20;
# Largeconf is used when the det space is very large but we have still no external file so we need to create one
largeconf = False
# Type of calculation
# TOTAL--> 1
# ELASTIC --> 2
# TOTAL ALIGNED --> 3
# ELASTIC ALIGNED --> 4
# TOTAL ELECTRON--> 5
# ELASTIC ELECTRON --> 6
# TOTAL J2 --> 7
# ELASTIC J2 --> 8
Type = 1
# Ouput name
mldfile = 'co_cas.mld'
punfile = 'co_cas.pun'
outfile = 'try_new_routine_2.dat'

readtwordm = False
file_read_twordm = 'QC-2RDM-AD.dat'
if not jeremyR and not hf and not readtwordm:

    # This routine reads the CIvectors and the configurations
    if molcas:
        logfile = 'molcas.log'
        civs, confs = mc.get_civs_and_confs(logfile, caspt2, extra_precision_molcas)
        print("civs")
        print(civs)
        print("confs")
        print(confs)

    elif bagel:
        print('will be implemented soon')
        confs = ['ab' * 18 + '00' * 2]
        civs = [1.000]

    elif molpro:
        civs, confs = td.twordmconst(closed, punfile)  # state1 and state2 should be used here

    Nmo_max = len(confs[0]) / 2
    print('Nmo_max:', Nmo_max)

elif not jeremyR and hf:
    civs = [1.000]
    # The number of occupied orbs or he configuration used must be specified by the user in hf
    norbs = 12
    confs = ['ab' * norbs]
    print(confs)
    Nmo_max = 12
else:
    # If the 2rdm or Civector is constructed in a bit-wise manner, the number or orbitals needs to be specified here by the user
    confs = 0
    civs = 0
    Nmo_max = 100
civs = np.array(civs)

if molcas:
    mldfile = 'molcas.rasscf.molden'
    gtos, atoms = mcmldreader.read_orbitals(mldfile, N=Nmo_max, decontract=True)
elif bagel:
    mldfile = 'orbitals.molden'
    Nmo_max = 20
    gtos, atoms = bgmldreader.read_orbitals(mldfile, N=Nmo_max, decontract=True)
elif molpro or readtwordm:
    Nmo_max = 100
    gtos, atoms, coeffs, mos, groupC,contr = mldreader.read_orbitals(mldfile, N=Nmo_max, decontract=False)
geom = atoms.geometry()
with open('options.dat', 'w') as f:
    f.write(str(np.size(atoms.atomic_numbers())) + '\n')
    for i in atoms.atomic_numbers():
        f.write(str(i) + ' ')
    f.write('\n')
    for i in range(np.size(atoms.atomic_numbers())):
        f.write(str(geom[i, :])[1:-1] + '\n')
    f.write(str(cutoffcentre) + '\n')
    f.write(str(cutoffz) + '\n')
    f.write(str(cutoffmd) + '\n')
    f.write(str(jeremyR) + '\n')
    f.write(str(mcci) + '\n')
    f.write(str(hf) + '\n')
    f.write(str(qmin) + ' ' + str(qmax) + ' ' + str(npoints) + ' \n')
    f.write(str(Type) + '\n')
    f.write(str(state1) + ' ' + str(state2) + '\n')
    f.write(outfile + '\n')
    f.write(str(molpro) + '\n')
    f.write(str(molcas) + '\n')
    f.write(str(bagel) + '\n')
    if molpro and not readtwordm and not hf:
        f.write(str(np.size(confs)) + '\n')

        f.write(str(len(confs[0])) + '\n')
        f.write(str(np.size(civs[0, :])) + '\n')

        if np.size(confs) < 1E5:
            bitwise = False
            print('Normal integration')
            for i in range(np.size(confs)):

                lst = wrap(confs[i].replace('a', '1').replace('b', '2'),
                           1)  # [*confs[i].replace('a', '1').replace('b', '2')]
                confs2 = np.asarray(lst, dtype=np.int64)

                g1 = [int(i) for i in lst]
                f.write(str(g1)[1:-1].replace(',', ' '))
                for j in range(np.size(civs[0, :])):
                    f.write(' ' + str(civs[i, j]) + ' ')
                f.write(str('\n'))
        else:
            print('bitwise integration')
            f.write('bitwise.dat')
            bitwise = True
            alpha = 0
            beta = 0
            with open('bitwise.dat', 'w') as m:
                for i in range(np.size(confs)):

                    lst = [*confs[i].replace('a', '1').replace('b', '2')]
                    confs2 = np.asarray(lst, dtype=np.int64)
                    alpha = 0
                    beta = 0
                    for j in range(np.size(confs2)):
                        if confs2[j] != 0:
                            if np.mod(j, 2) == 0:
                                alpha = alpha + 2 ** (j / 2)
                            else:
                                beta = beta + 2 ** ((j - 1) / 2)
                    m.write(str(i) + ' ' + str(civs[i, state1 - 1]) + ' ' + str(civs[i, state2 - 1]) + ' ' + str(
                        int(alpha)) + ' ' + str(int(beta)) + '\n')
        f.write(str(bitwise) + '\n')
    elif bagel and not hf:
        #  f.write('bagel')
        f.write(str(np.size(confs)) + '\n')
        print(str(np.size(civs[:, 0])))
        f.write(str(len(confs[0])) + '\n')
        f.write(str(np.size(civs[0, :])) + '\n')
        lst = [*confs[0].replace('a', '1').replace('b', '2')]
        confs2 = np.asarray(lst, dtype=np.int64)
        f.write(str(confs2)[1:-1])
        f.write(' ' + str(1.00000) + ' ')
        f.write(str('\n'))
    elif bagel and hf:
        f.write(str(1) + '\n')
        f.write(str(40) + '\n')
        f.write(str(1) + '\n')
        lst = [*confs[0].replace('a', '1').replace('b', '2')]
        confs2 = np.asarray(lst, dtype=np.int64)
        f.write(str(confs2)[1:-1])
        f.write(' ' + str(1.00000) + ' ')
        f.write(str('\n'))

    elif molcas and not hf:
        f.write(str(np.size(confs)) + '\n')
        print(str(np.size(civs[:, 0])))
        f.write(str(len(confs[0])) + '\n')
        f.write(str(np.size(civs[0, :])) + '\n')
        #  f.write('molcas')
        for i in range(np.size(confs)):
            lst = wrap(confs[i].replace('a', '1').replace('b', '2'), 1)
            confs2 = [int(i) for i in lst]
            f.write(str(confs2)[1:-1].replace(',', ''))
            for j in range(np.size(civs[0, :])):
                f.write(' ' + str(civs[i, j]) + ' ')
            f.write(str('\n'))
    elif readtwordm:
        f.write('readtwordm' + '\n')
        f.write(file_read_twordm)
    elif hf:
        f.write(str(np.size(confs)) + '\n')

        f.write(str(len(confs[0])) + '\n')

        f.write(str(1) + '\n')

        lst = [*confs[0].replace('a', '1').replace('b', '2')]
        confs2 = np.asarray(lst, dtype=np.int64)
        f.write(str(confs2)[1:-1])
        f.write(' ' + str(1.00000) + ' ')
        f.write(str('\n'))

Nmo_max = 100

print(geom)
xx = gtos.x
yy = gtos.y
zz = gtos.z
l = gtos.l
m = gtos.m
n = gtos.n
ga = gtos.ga
group = gtos.group

mmod = np.transpose(gtos.mo)

l = np.asarray(l)
m = np.asarray(m)
n = np.asarray(n)
ga = np.asarray(ga)
mmod = np.asarray(mmod, dtype=np.float64)

print(np.argsort(group))



with open('basis.dat', 'w') as f:
    f.write(str(np.size(l)) + '\n')
    for i in range(np.size(l)):
        f.write(str(xx[i]) + ' ' + str(yy[i]) + ' ' + str(zz[i]) + ' ' + str(ga[i]) + ' ' + str(l[i]) + ' ' + str(
            m[i]) + ' ' + str(n[i]) + ' ' + str(group[i]) + ' '+ str(contr[i])+'\n')
with open('MOs.dat', 'w') as f:
    f.write(str(np.size(mmod[:, 0])) + ' ' + str(np.size(mmod[0, :])) + '\n')
    for i in range(np.size(mmod[:, 0])):
        for j in range(np.size(mmod[0, :])):
            f.write(str(mmod[i, j]) + ' ')
        f.write('\n')
with open('coeffs.dat', 'w') as f:
    f.write(str(np.size(l)) + '\n')
    for i in range(np.size(l)):
        f.write(str(coeffs[i]) + '\n')
    f.write(str(np.size(groupC)) + '\n')
    count = 1
    for i in range(np.size(groupC)):
        f.write(str(count) + ' ' + str(count + groupC[i] - 1) + ' ' + str(groupC[i]) + '\n')
        count = count + groupC[i]
