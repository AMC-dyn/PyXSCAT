import numpy as np
import molproinp_out as mp
import ReorderGTOS as Rdgt
import twordm_red as td
import twordm as td2
from factd import factd
from pzerotable import pzerotablenew as pzerot
# from MDcoeffs import md_table_gen as mtg
# from integral_ijkr import integral_k_ijkr as intk
# from integral_ijkr_vect import integral_k_ijkr as intk_vec
# from integral_ijkr_pzero import integral_k_ijkr_0_case as intkzero
from unsortunique import uunique
from integrals_wrapper import integrals_ijkr
import time
import molden_reader_nikola as mldreader


# Rotational averaged  caclulation of total scattering

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
    x, y, z = mp.create_input()
    inputime2 = time.time()
    # This routine creates the molpro input, a punch file with the CI vectors and a molden file where
    # the basis and MO coefficients
    # are read
   # ga, ci, realnum, m, l, n, mos, monums, actives, total, angpart, xx, yy, zz = mp.outputreading(x, y, z)
   # fulltable = np.array(np.vstack((ga, l, m, n, xx, yy, zz)))  # Note: is this table correct?
   # dummy, ipos, irep = uunique(fulltable)

    # norma = np.power(np.divide(2, np.pi), 0.75)
    # normb = np.zeros(l.size)
    # normc = np.zeros(l.size)
    # normd = np.zeros(l.size)
    # for i in range(l.size):
    #     normb[i] = np.power(2, l[i] + m[i] + n[i])
    #     normc[i] = np.power(ga[i],
    #                         np.divide(np.multiply(2, l[i]) + np.multiply(2, m[i]) + np.multiply(2, n[i]) + 3, 4.0))
    #     normd[i] = np.multiply(np.multiply(factd(np.multiply(2, l[i]) - 1), factd(np.multiply(2, m[i]) - 1)),
    #                            factd(np.multiply(2, n[i]) - 1))
    # nrm = norma * normb * normc / normd ** 0.5
    print('input time', inputime2 - inputime1, 's')

    mldfile = 'molpro.mld'
    Nmo_max = 9

    gtos = mldreader.read_orbitals(mldfile, N=Nmo_max, decontract=True)

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

    # mos = np.reshape(mos, (total, max(realnum)))
    # mos = mos[:actives, :]
    # realnum = np.asarray(realnum)
    # mos = mos[:, realnum - 1]



    # The original GTOs are reordered and MOS is changed,converted in a matrix and reordered
    # l, m, n, ga, ci, mos, angpart = Rdgt.reorder(l, m, n, ga, ci, mos, angpart)

    print("Cutoff values are specified by default as 0.01, 1E-9, 1E-20\n")
    # condit = input("Do you want to continue Y/N?")
    x = True
    if x:
        # cut off epsilon; if H < epsilon, use P0 cases
        cutoffcentre = 0.01  # suggested: cutoffcentre = 0.01;
        # the cutoff for the Z integral
        cutoffz = 1e-9  # suggested: cutoffz = 1E-9;
        # the cutoff for the product of the MD coefficients
        cutoffmd = 1e-30  # suggested: cutoffmd = 1E-20;
    else:
        cutoffcentre = input("Input cut off epsilon; if H < epsilon, use P0 cases")
        cutoffmd = input("Input the cutoff for the product of the MD coefficients")
        cutoffz = input("Input the cutoff for the Z integral")
    # reshape q to a column vector
    q = np.linspace(0.000001, 3, 2)
    # set q t0 E-10 if q = 0
    if q[0] < 1E-10:
        q[0] = 1E-10

    # reading the 2-particle RDM from MOLPRO output

    civs, confs = td.twordmconst()  # state1 and state2 should be used here
    # dumm2,dumm3 = td2.twordmconst()
    print('twordm constructed')
    # normalisation of the GTOs

    # combine the normalisation, MO, and contraction coefficients

   # mmod = mos * ci * nrm

    # total number of primitive GTOs
    ncap = l.size
    print('The number of primitive GTOs before reduction: ', str(ncap))

    # generate table with GTO descriptors

    print(l)
    # dummy,ipos,irep = unique(full_table,'rows','stable')

    # ipos = np.sort(ipos)
    # irep = np.sort(irep)
    # print(irep)
    # print(dummy)
    # print(np.size(irep), ncap, ipos.size, "no pasaran, las dimensiones, no pasaran")
    # reduce the matrices to include only the unique cGTOs
    # ga = ga[ipos]
    # l = np.asarray(l[ipos], dtype=np.int32)  # AMBIGUOUS VARIABLE NAME?
    # # print(ipos, irep, "fascistas")
    # m = m[ipos]
    # n = n[ipos]
    # xx = xx[ipos]
    # yy = yy[ipos]
    # zz = zz[ipos]
    # angpart = angpart[ipos]
    #
    # # total number of primitive GTOs
    # ncap = l.size

    print('The number of primitive GTOs after reduction: ', str())

    # print("p0 first point:", p0matrix[:, 1, 0, 0])
    # m1 = np.asarray(mat[:, 0], dtype=np.int64) + 1
    # m2 = np.asarray(mat[:, 1], dtype=np.int64) + 1
    # m3 = np.asarray(mat[:, 2], dtype=np.int64) + 1
    # m4 = np.asarray(mat[:, 3], dtype=np.int64) + 1

    potime1 = time.time()
    # p0matrix = np.zeros(
    #     (np.size(q), 4 * np.max(l) + 2, 4 * np.max(l) + 2, 4 * np.max(l) + 2), dtype=np.float64)
    # for i in range(4 * np.max(l) + 2):
    #     for j in range(4 * np.max(l) + 2):
    #         for k in range(4 * np.max(l) + 2):
    #             p0matrix[:, i, j, k] = pzerot(i, j, k, q)
    potime2 = time.time()
    print('p0time ', potime2 - potime1, 's')
    # Calculation of the prexponential factors Z and Z2 for all N^2 GTO products and N_orb
    Ngto = np.size(l)
    ng = max(group)

    # l1 = np.zeros((Ngto, len(irep)))
    # ln1 = np.zeros(Ngto)
    #
    # for i in range(Ngto):
    #     iduplicates = 1 + np.argwhere(irep == irep[ipos[i]])
    #
    #     iduplicates = np.asarray(iduplicates).flatten()
    #
    #     l1[i, :np.size(iduplicates)] = iduplicates
    #
    #     ln1[i] = np.size(iduplicates)
    #
    # # Dealing with the orbitals of a similar type together
    # dummy, apos, arep = np.unique(angpart, return_index=True, return_inverse=True)  # Note: is "axis=0" correct?
    #
    # nnew = apos.size
    # print('The number of GTOs after compression is: ', str(nnew))
    #
    # listofdups2 = np.zeros((nnew, len(arep)))
    # listofnumbers2 = np.zeros(nnew)
    #
    # for i in range(nnew):
    #     iduplicates = 1 + np.argwhere(arep == arep[apos[i]])
    #
    #     iduplicates = np.asarray(iduplicates).flatten()
    #
    #     listofdups2[i, :np.size(iduplicates)] = iduplicates
    #
    #     listofnumbers2[i] = np.size(iduplicates)
    #
    # # definition of the new size
    # ncap = nnew
    #
    # # Loops over 4xpGTOs
    # # Strictly non-diagonal elements with all GTOs being different
    # # int_res = np.zeros(q.size)
    #
    # ipos = np.array(ipos, dtype=np.int) + 1
    #
    # apos = np.array(apos, dtype=np.int) + 1
    #
    maxl = max(l)
    nq = np.size(q)
    #
    # nipos = np.size(ipos)
    #
    # napos = ncap
    confs2 = np.zeros((np.size(confs), len(confs[0])), dtype=np.int64)
    counts = 0
    for i in confs:
        lst = []
        lst = [*i.replace('a', '1').replace('b', '2')]
        confs2[counts, :] = np.asarray(lst, dtype=np.int64)
        counts = counts + 1

    # print(np.max(p0matrix))
    # px, py, pz, dx, dy, dz, z1, z2, e12, ll = integrals_ijkr.variables_total_3(maxl, ipos, nnn2, apos, nnew, ga,
    #                                                                          l, m, n,
    #                                                                          xx, yy, zz, mmod, m1, m2, m3, m4,
    #                                                                          nnn3,
    #                                                                          total, q, nq, listofdups1,
    #                                                                          listofnumbers1,
    #                                                                          listofdups2, listofnumbers2)
    # print(px[1, 1])
    # resultado2 = np.zeros(nq)
    # resultado2 = integrals_ijkr.integration(ncap, px, py, pz, ll, p0matrix, dx, dy, dz, z1, z2, apos, cutoffz,
    #                                         cutoffmd,
    #                                         cutoffcentre, q, e12)
    tic2 = time.time()
    print(np.size(group))
    print('time for readers in python', tic2 - tic1, 's')
    resultado2 = integrals_ijkr.total_scattering_calculation(maxl, Ngto, ng,
                                                             ga, l, m, n, xx, yy, zz, mmod,
                                                             q, nq,
                                                             group,
                                                             cutoffz, cutoffmd, cutoffcentre, confs2, civs)
    print(resultado2)

    return 1


tic2 = time.time()

re = main()
toc = time.time()
print(toc - tic2)
# print('whats the error')
