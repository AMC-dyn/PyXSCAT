import numpy as np
import molproinp_out as mp
import twordm_red as td
import twordm as td2
from integrals_wrapper import main_calculation_mod as main_calculation
import time
import molden_reader_nikola as mldreader
import Molcas_output_reader as Molcas
import scipy.io as sci
from PyQt5.QtCore import *
from PyQt5.QtGui import *
from PyQt5.QtWidgets import *
import sys
import pandas as pd


# Rotational averaged  caclulation of total scattering
def integertobinary(N):
    binary = []
    for i in range(0, 32):
        k = int(np.mod(N, 2))
        binary.append(k)
        N = N / 2

    return binary


def read_fci(file):
    CIV = False
    alphas = []
    betas = []
    civectors = []

    with open(file, 'r') as f:

        for lines in f:
            if 'CI Vector' in lines:
                CIV = True
                print('CIV FOUND')

            if CIV and '/EOF' not in lines and 'CI Vector' not in lines:

                N = lines.strip().split()
                vec = np.array(N[1:], dtype=int)

                civ = np.array(N[0], dtype=np.float64)
                alpha, beta = orbitals_to_integers(vec)

                civectors.append(civ)
                alphas.append(alpha)
                betas.append(beta)
            elif '/EOF' in lines:
                CIV = False
                break
    civectors = np.asarray(civectors)
    alphas = np.asarray(alphas)
    betas = np.asarray(betas)

    return civectors, alphas, betas


def orbitals_to_integers(vec):
    alpha = 0
    beta = 0
    comp = 0
    a = True
    b = False
    for i in vec:
        if i < comp:
            b = True
            a = False
        if a:
            alpha += 2 ** (i - 1)

        if b:
            beta += 2 ** (i - 1)

        comp = i
    return alpha, beta


def main():
    # This function creates the input and produces the molpro output or just does not run if it is specified in abinitio.dat
    # Also it copies the outputs into molpro.mld and molpro.pun automatically

    mp.create_input()
    mldfile = 'molpro.mld'

    # If there is an external file with CIvecs or 2rdms JeremyR==True
    jeremyR = True
    fileJeremy = 'configurations_bit.dat'
    # If we convert form a MCCI calculation to a bitwise operation mcci==True
    mcci = False
    # If we have a HF calculation and our Civector is represented by a single determinant hf==True
    hf = False
    # States involved
    state1 = 1
    state2 = 2
    # Largeconf is used when the det space is very large but we have still no external file so we need to create one
    largeconf = True
    #Type of calculation
    #TOTAL--> 1
    #ELASTIC --> 2
    #TOTAL ALIGNED --> 3
    #ELASTIC ALIGNED --> 4
    # TOTAL ELECTRON--> 5
    # ELASTIC ELECTRON --> 6
    # TOTAL J2 --> 7
    # ELASTIC J2 --> 8
    Type=9
    # Ouput name (is a mat file, needs to be changed to make it general)
    nameoffile = 'NO_CMS_12_j1.mat'

    if not jeremyR and not hf:

        # This routine reads the CIvectors and the configurations
        civs, confs = td.twordmconst()  # state1 and state2 should be used here

        # confs, civs = Molcas.reading_molcas("big_AS/molcas.log", 1)

        Nmo_max = len(confs[0]) / 2
        print('Nmo_max:', Nmo_max)

    elif not jeremyR and hf:
        civs = 1.000
        # The number of occupied orbs or he configuration used must be specified by the user in hf
        norbs = 7
        confs = ['ab' * norbs]
        print(confs)
        Nmo_max = 18
    else:
        # If the 2rdm or Civector is constructed in a bit-wise manner, the number or orbitals needs to be specified here by the user
        Nmo_max = 14

    print('Max_nmos,', Nmo_max)
    # gtos, atoms = mldreader.read_orbitals('2andres/molcas.rasscf.molden.1', N=Nmo_max, decontract=True)

    gtos, atoms = mldreader.read_orbitals(mldfile, N=Nmo_max, decontract=True)
    geom = atoms.geometry()
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

    # Reading and transforming the molden part (basis, geometry and MOs)

    print("Cutoff values are specified by default as 0.01, 1E-9, 1E-20\n")
    # condit = input("Do you want to continue Y/N?")
    x = True
    if x:
        # cut off epsilon; if H < epsilon, use P0 cases
        cutoffcentre = 0.01  # suggested: cutoffcentre = 0.01;
        # the cutoff for the Z integral
        cutoffz = 1e-15  # suggested: cutoffz = 1E-9;
        # the cutoff for the product of the MD coefficients
        cutoffmd = 1e-20  # suggested: cutoffmd = 1E-20;
    else:
        cutoffcentre = input("Input cut off epsilon; if H < epsilon, use P0 cases")
        cutoffmd = input("Input the cutoff for the product of the MD coefficients")
        cutoffz = input("Input the cutoff for the Z integral")
    # reshape q to a column vector
    q = np.linspace(0.00000000001, 10, 200)
    # set q t0 E-10 if q = 0
    if q[0] < 1E-10:
        q[0] = 1E-10

    # total number of primitive GTOs
    ncap = l.size
    print('The number of primitive GTOs before reduction: ', str(ncap))

    # generate table with GTO descriptors

    # Calculation of the prexponential factors Z and Z2 for all N^2 GTO products and N_orb
    Ngto = np.size(l)
    ng = max(group)
    maxl = max(l)
    nq = np.size(q)

    ###### FIRST CASE###################
    if largeconf and not jeremyR:
        jeremyR = True
        fileJeremy = 'configurations_bit.dat'
        confs2 = np.zeros((np.size(confs), len(confs[0])), dtype=np.int64)

        counts = 0
        for i in confs:
            lst = []
            lst = [*i.replace('a', '1').replace('b', '2')]
            confs2[counts, :] = np.asarray(lst, dtype=np.int64)
            counts = counts + 1

        with open('configurations_bit.dat', 'w') as f:
            for i in range(counts - 1):
                var1 = confs2[i, :]
                alpha = 0
                beta = 0

                for j in range(len(var1)):
                    if var1[j] != 0:
                        if np.mod(j, 2) == 0:
                            alpha = alpha + 2 ** (j / 2)
                        else:
                            beta = beta + 2 ** ((j - 1) / 2)




    ## There is not file specified and
    if not jeremyR:
        confs2 = np.zeros((np.size(confs), len(confs[0])), dtype=np.int64)
        counts = 0
        for i in confs:
            lst = []
            lst = [*i.replace('a', '1').replace('b', '2')]
            confs2[counts, :] = np.asarray(lst, dtype=np.int64)
            counts = counts + 1

        q_alig, resultado2 = main_calculation.total_scattering_calculation(Type, atoms.atomic_numbers(), geom, state1, state2, maxl,
                                                                           Ngto, ng,
                                                                           ga, l, m, n, xx, yy, zz, mmod,
                                                                           q, nq,
                                                                           group,
                                                                           cutoffz, cutoffmd, cutoffcentre, confs2,
                                                                           civs)

    elif jeremyR:
        count = 0
        fci = False
        read2rdm = False
        if not fileJeremy:
            fileJeremy = 'configurations_bit.dat'

        if fci:
            civs, alphas, betas = read_fci(fileJeremy)
            test = np.sum(civs ** 2)

            # matrix2 = np.asarray([alphas, betas])
            # matrix2 = np.sort(matrix2, axis=0)
            # print(matrix2[:, 4763253])

            index = np.where(abs(civs) >= 1E-8)
            print('the number of elements selected is ', np.size(index))
            alphas = alphas[index]
            betas = betas[index]
            civs = civs[index]
            newdat, ipos1, irep1 = np.unique(alphas, return_index=True, return_inverse=True)
            newdat2, ipos2, irep2 = np.unique(betas, return_index=True, return_inverse=True)

            nnn = np.ones(np.size(newdat), dtype=np.int)
            start = np.ones(np.size(newdat), dtype=np.int)
            end = np.ones(np.size(newdat), dtype=np.int)
            nnn2 = np.ones(np.size(newdat), dtype=np.int)
            start_2 = np.ones(np.size(newdat2), dtype=np.int)
            end_2 = np.ones(np.size(newdat2), dtype=np.int)
            for i in range(np.size(newdat)):
                var = np.where(alphas == newdat[i])[0]
                nnn[i] = np.size(var)
                start[i] = int(var[0] + 1)
                end[i] = int(var[-1] + 1)
            for i in range(np.size(newdat2)):
                var = np.where(betas == newdat2[i])[0]
                nnn2[i] = np.size(var)
                start_2[i] = int(var[0] + 1)
                end_2[i] = int(var[-1] + 1)

            irep1 = nnn

            if np.size(newdat) < np.size(newdat2):
                newarray = np.zeros(np.size(newdat2) - np.size(newdat))
                newdat = np.append(newdat, newarray)
            elif np.size(newdat) > np.size(newdat2):
                newarray = np.zeros(np.size(newdat) - np.size(newdat2))
                newdat2 = np.append(newdat2, newarray)
            newdat = np.transpose(np.asarray([newdat, newdat2]))

            test = np.sum(civs[abs(civs) != 1E-6] ** 2)

            print('dataframe created')
            f = open('fci_calc.dat', 'w')

            for i in range(0, len(civs)):
                if abs(civs[i]) >= 1E-8:
                    f.write(str(i) + ' ' + str(civs[i]) + ' ' + str(alphas[i]) + ' ' + str(betas[i]) + '\n')
                    count = count + 1
            f.close()

            fileJeremy = 'fci_calc.dat'
            # count = len(civs)
            ordering1 = np.zeros(np.size(alphas))
            ordering2 = np.zeros(np.size(alphas))
            q_alig, resultado2 = main_calculation.total_scattering_calculation_2(type, atoms.atomic_numbers(), geom, state1,
                                                                                 state2,
                                                                                 maxl,
                                                                                 Ngto, ng,
                                                                                 ga, l, m, n, xx, yy, zz, mmod,
                                                                                 q, nq,
                                                                                 group,
                                                                                 cutoffz, cutoffmd, cutoffcentre,
                                                                                 fileJeremy, count, newdat, start,
                                                                                 end, start_2, end_2, fci, nnn,
                                                                                 nnn2, ordering1, ordering2, read2rdm,
                                                                                 mcci)

        else:
            if not read2rdm:
                confbin1 = []
                confbin2 = []
                civs = []
                f = open(fileJeremy, 'r')
                for line in f:
                    NN = line.strip().split()
                    confbin1.append(int(NN[3]))
                    confbin2.append(int(NN[4]))
                    count += 1
                f.close()
                if count <= 300:
                    confbin1 = []
                    confbin2 = []
                    f = open(fileJeremy, 'r')
                    for line in f:
                        NN = line.strip().split()

                        civs.append(float(NN[1]))
                        confbin1.append(integertobinary(int(NN[2])))
                        confbin2.append(integertobinary(int(NN[3])))

                    confs2 = np.zeros((count, len(confbin1[0] * 2)))
                    confbin1 = np.asarray(confbin1)
                    confbin2 = np.asarray(confbin2)
                    confbin2 = confbin2 * 2
                    civs = np.asarray(civs, dtype=np.float64)
                    print('norm of the CI', sum(civs ** 2))
                    civs = civs / np.sqrt(sum(civs ** 2))
                    for i in range(0, count):
                        for j in range(0, 32):
                            confs2[i, 2 * j] = confbin1[i, j]
                            confs2[i, 2 * j + 1] = confbin2[i, j]
                    print('confs2', len(confs2))

                    q_alig, resultado2 = main_calculation.total_scattering_calculation(Type, atoms.atomic_numbers(), geom,
                                                                                       state1,
                                                                                       state2,
                                                                                       maxl,
                                                                                       Ngto, ng,
                                                                                       ga, l, m, n, xx, yy, zz, mmod,
                                                                                       q, nq,
                                                                                       group,
                                                                                       cutoffz, cutoffmd, cutoffcentre,
                                                                                       confs2,
                                                                                       civs)
                else:

                    ordering1 = 0
                    ordering2 = 0
                    # print(np.asarray(confbin1)[ordering1])
                    newdat, ipos, irep1 = np.unique(confbin1, return_index=True, return_inverse=True)
                    newdat2, ipos, irep2 = np.unique(confbin2, return_index=True, return_inverse=True)
                    nnn = np.ones(np.size(newdat), dtype=np.int)
                    start = np.ones(np.size(newdat), dtype=np.int)
                    end = np.ones(np.size(newdat), dtype=np.int)
                    nnn2 = np.ones(np.size(newdat2), dtype=np.int)
                    start_2 = np.ones(np.size(newdat2), dtype=np.int)
                    end_2 = np.ones(np.size(newdat2), dtype=np.int)

                    # for i in range(np.size(newdat)):
                    #     var = np.where(np.asarray(confbin1)[ordering1] == newdat[i])[0]
                    #     nnn[i] = np.size(var)
                    #     start[i] = int(var[0] + 1)
                    #     end[i] = int(var[-1] + 1)
                    # for i in range(np.size(newdat2)):
                    #     var = np.where(np.asarray(confbin2)[ordering2] == newdat2[i])[0]
                    #
                    #     nnn2[i] = np.size(var)
                    #     start_2[i] = int(var[0] + 1)
                    #     end_2[i] = int(var[-1] + 1)
                    print(np.size(newdat), np.size(newdat2))
                    if np.size(newdat) < np.size(newdat2):
                        newarray = np.zeros(np.size(newdat2) - np.size(newdat))
                        newdat = np.append(newdat, newarray)
                    elif np.size(newdat) > np.size(newdat2):
                        newarray = np.zeros(np.size(newdat) - np.size(newdat2))
                        newdat2 = np.append(newdat2, newarray)

                    print(np.size(newdat), np.size(newdat2))
                    newdat = np.transpose(np.asarray([newdat, newdat2]))

                    print(np.shape(newdat))

                    q_alig, resultado2 = main_calculation.total_scattering_calculation_2(Type, atoms.atomic_numbers(),
                                                                                         geom, state1,
                                                                                         state2,
                                                                                         maxl,
                                                                                         Ngto, ng,
                                                                                         ga, l, m, n, xx, yy, zz, mmod,
                                                                                         q, nq,
                                                                                         group,
                                                                                         cutoffz, cutoffmd,
                                                                                         cutoffcentre,
                                                                                         fileJeremy, count, newdat,
                                                                                         start,
                                                                                         end, start_2, end_2, fci, nnn,
                                                                                         nnn2, ordering1, ordering2,
                                                                                         read2rdm, mcci)

            else:
                newdat = 0
                start = 0
                end = 0
                start_2 = 0
                count = 0
                fci = False
                nnn = 0
                nnn2 = 0
                ordering1 = 0
                ordering2 = 0
                end_2 = 0
                with open(fileJeremy, 'r') as f:
                    for lines in f:
                        count = count + 1
                q_alig, resultado2 = main_calculation.total_scattering_calculation_2(Type, atoms.atomic_numbers(), geom, state1,
                                                                                     state2,
                                                                                     maxl,
                                                                                     Ngto, ng,
                                                                                     ga, l, m, n, xx, yy, zz, mmod,
                                                                                     q, nq,
                                                                                     group,
                                                                                     cutoffz, cutoffmd, cutoffcentre,
                                                                                     fileJeremy, count, newdat, start,
                                                                                     end, start_2, end_2, fci, nnn,
                                                                                     nnn2, ordering1, ordering2,
                                                                                     read2rdm, mcci)

    # tic2 = time.time()
    # print(np.size(group))
    # print('Angular momenta red', ng)
    # print('time for readers in python', tic2 - tic1, 's')

    print(resultado2)
    # print(q)
    return resultado2, q, q_alig,nameoffile


tic2 = time.time()

res, q, q_alig,nameoffile = main()
toc = time.time()

sci.savemat(nameoffile, {'q': q, 'I': res})

print(nameoffile)
print(toc - tic2)
