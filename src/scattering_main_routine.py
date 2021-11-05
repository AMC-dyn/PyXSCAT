import numpy as np
import molproinp_out as mp
import twordm_red as td
import twordm as td2
from integrals_wrapper import main_calculation_mod as main_calculation
import time
import molden_reader_nikola as mldreader
import scipy.io as sci


# Rotational averaged  caclulation of total scattering
def integertobinary(N):
    binary = []
    for i in range(0, 32):
        k = int(np.mod(N, 2))
        binary.append(k)
        N = N / 2

    return binary


def main():
    tic1 = time.time()
    # INPUT:
    # mldfile       string          input molden file name
    # outifle       string          molpro output file name
    # state1        int             state1 to consider (NOT USED!)
    # state2        int             if state2 == state1, total; else coherences (NOT USED!)
    # nstates       int             the number of states in the state-averaged calculation found in the molpro file
    #                               (NOT USED!)
    # q[0:nq]       vec             scattering vector in ATOMIC UNITS; q[0] must be different from 0
    # cutoff[0:2]   vec             the set of cutoff values: centre, Zijkr integral, MD part of Z

    # OUTPUT:
    # tsi[0:nq]     column vec      scattering signal

    # We need to make the molpro execution voluntary, otherwise the punch file and the molden file
    # are copied into molpro.pun and molpro.mld
    inputime1 = time.time()
    mp.create_input()
    inputime2 = time.time()

    print('input time', inputime2 - inputime1, 's')

    mldfile = 'molpro.mld'
    Nmo_max = 32
    jeremyR = True
    if not jeremyR:
        civs, confs = td.twordmconst()  # state1 and state2 should be used here
        Nmo_max = len(confs[:][0]) / 2
    else:
        Nmo_max = 32
    print('Max_nmos,', Nmo_max)
    gtos, atoms = mldreader.read_orbitals(mldfile, N=Nmo_max, decontract=True)

    geom = atoms.geometry()
    print(atoms.atomic_numbers())

    print(np.size(gtos.ga))

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

    print("Cutoff values are specified by default as 0.01, 1E-9, 1E-20\n")
    # condit = input("Do you want to continue Y/N?")
    x = True
    if x:
        # cut off epsilon; if H < epsilon, use P0 cases
        cutoffcentre = 0.01  # suggested: cutoffcentre = 0.01;
        # the cutoff for the Z integral
        cutoffz = 1e-30  # suggested: cutoffz = 1E-9;
        # the cutoff for the product of the MD coefficients
        cutoffmd = 1e-30  # suggested: cutoffmd = 1E-20;
    else:
        cutoffcentre = input("Input cut off epsilon; if H < epsilon, use P0 cases")
        cutoffmd = input("Input the cutoff for the product of the MD coefficients")
        cutoffz = input("Input the cutoff for the Z integral")
    # reshape q to a column vector
    q = np.linspace(0.00000000001, 10, 100)
    # set q t0 E-10 if q = 0
    if q[0] < 1E-10:
        q[0] = 1E-10

    # reading the 2-particle RDM from MOLPRO output

    # mattry, total = td2.twordmconst()

    print('twordm constructed')

    # total number of primitive GTOs
    ncap = l.size
    print('The number of primitive GTOs before reduction: ', str(ncap))

    # generate table with GTO descriptors

    print(l)

    potime1 = time.time()

    potime2 = time.time()
    print('p0time ', potime2 - potime1, 's')
    # Calculation of the prexponential factors Z and Z2 for all N^2 GTO products and N_orb
    Ngto = np.size(l)
    ng = max(group)
    maxl = max(l)
    nq = np.size(q)

    if not jeremyR:
        confs2 = np.zeros((np.size(confs), len(confs[0])), dtype=np.int64)
        counts = 0
        for i in confs:
            lst = []
            lst = [*i.replace('a', '1').replace('b', '2')]
            confs2[counts, :] = np.asarray(lst, dtype=np.int64)
            counts = counts + 1

        q_alig, resultado2 = main_calculation.total_scattering_calculation(2, atoms.atomic_numbers(), geom, 1, 1, maxl,
                                                                           Ngto, ng,
                                                                           ga, l, m, n, xx, yy, zz, mmod,
                                                                           q, nq,
                                                                           group,
                                                                           cutoffz, cutoffmd, cutoffcentre, confs2,
                                                                           civs)






    elif jeremyR:
        count = 0
        fileJeremy = 'O3/civ_out_0.0001'
        confbin1 = []
        confbin2 = []
        civs = []
        f = open(fileJeremy, 'r')
        for line in f:

            count += 1
        f.close()
        if count <= 10:
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

            q_alig, resultado2 = main_calculation.total_scattering_calculation(1, atoms.atomic_numbers(), geom, 1, 1,
                                                                               maxl,
                                                                               Ngto, ng,
                                                                               ga, l, m, n, xx, yy, zz, mmod,
                                                                               q, nq,
                                                                               group,
                                                                               cutoffz, cutoffmd, cutoffcentre, confs2,
                                                                               civs)
        else:

            q_alig, resultado2 = main_calculation.total_scattering_calculation_2(1, atoms.atomic_numbers(), geom, 1, 1,
                                                                               maxl,
                                                                               Ngto, ng,
                                                                               ga, l, m, n, xx, yy, zz, mmod,
                                                                               q, nq,
                                                                               group,
                                                                               cutoffz, cutoffmd, cutoffcentre, fileJeremy,count)
    #tic2 = time.time()
    #print(np.size(group))
    #print('Angular momenta red', ng)
    #print('time for readers in python', tic2 - tic1, 's')

    print(resultado2)
    #print(q)
    return resultado2, q, q_alig


tic2 = time.time()

res, q, q_alig = main()
toc = time.time()
sci.savemat('CO_0p0005.mat', {'q': q, 'I': res})
print(toc - tic2)
