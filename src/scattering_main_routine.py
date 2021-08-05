import numpy as np
import molproinp_out as mp
import ReorderGTOS as rdgt
import twordm as td
from factd import factd
from pzerotable import pzerotablenew as pzerot
from MDcoeffs import md_table_gen as mtg
from integral_ijkr import integral_k_ijkr as intk
from integral_ijkr_vect import integral_k_ijkr as intk_vec
from integral_ijkr_pzero import integral_k_ijkr_0_case as intkzero
from unsortunique import uunique
from integrals_wrapper import integrals_ijkr


# Rotational averaged  caclulation of total scattering

def main():
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
    x, y, z = mp.create_input()

    x = [0.0, 0.0]
    y = [0.0, 0.0]
    z = [0.0, 1.09]
    # This routine creates the molpro input, a punch file with the CI vectors and a molden file where the basis and MO coefficients
    # are read
    ga, ci, realnum, m, l, n, mos, monums, actives, total, angpart, xx, yy, zz = mp.outputreading(x, y, z)

    # mos = np.reshape(mos, (total, max(realnum)))
    # mos = mos[:actives, :]
    # realnum = np.asarray(realnum)
    # mos = mos[:, realnum - 1]

    # Change everything to atomic units (originally in Angstroms)
    xx = xx / 0.529177210920
    yy = yy / 0.529177210920
    zz = zz / 0.529177210920

    # The original GTOs are reordered and MOS is changed,converted in a matrix and reordered
    l, m, n, ga, ci, mos, angpart = rdgt.reorder(l, m, n, ga, ci, mos, angpart)

    # mldfile = sys.argv[1]
    # outfile = sys.argv[2]
    # state1 = sys.argv[3]  # Note: NOT USED ANYWHERE!
    # state2 = sys.argv[4]  # Note: NOT USED ANYWHERE!
    # nstates = sys.argv[5]  # Note: NOT USED ANYWHERE!
    # q = float(sys.argv[6])

    print("Cutoff values are specified by default as 0.01, 1E-9, 1E-20\n")
    # condit = input("Do you want to continue Y/N?")
    x = True
    if x:
        # cut off epsilon; if H < epsilon, use P0 cases
        cutoffcentre = 0.01  # suggested: cutoffcentre = 0.01;
        # the cutoff for the Z integral
        cutoffz = 1e-9  # suggested: cutoffz = 1E-9;
        # the cutoff for the product of the MD coefficients
        cutoffmd = 1e-20  # suggested: cutoffmd = 1E-20;
    else:
        cutoffcentre = input("Input cut off epsilon; if H < epsilon, use P0 cases")
        cutoffmd = input("Input the cutoff for the product of the MD coefficients")
        cutoffz = input("Input the cutoff for the Z integral")
    # hola
    # def tot_rot_av(mldfile,outfile,state1,state2,nstates,q,cutoffs):

    # # cut off epsilon; if H < epsilon, use P0 cases
    # cutoffcentre = cutoffs[0]  # suggested: cutoffcentre = 0.01;
    # # the cutoff for the Z integral
    # cutoffz = cutoffs[1]  # suggested: cutoffz = 1E-9;
    # # the cutoff for the product of the MD coefficients
    # cutoffmd = cutoffs[2]  # suggested: cutoffmd = 1E-20;

    # reshape q to a column vector
    q = np.linspace(0.000001, 3, 100)
    # q = np.reshape(q, (q.size))
    # q = np.reshape(q, (q.size, 1))

    # set q t0 E-10 if q = 0
    if q[0] < 1E-10:
        q[0] = 1E-10

    # reading the 2-particle RDM from MOLPRO output
    mat, total = td.twordmconst()  # state1 and state2 should be used here

    nmocutoff = np.max(np.max(mat))
    # atoms, ga, xx, yy, zz, l, m, n, c, mos, dummy, angpart, ntc, reps, matoms = mldread(mldfile, nmocutoff)
    # Note: CHECK AND IMPORT FUNCTION "mldread"; CHANGE CHARACTER "M"

    # reording the GTOs uising the angular momenta
    # l, m, n, c, mos, ga, angpart = reorder(l, m, n, c, mos, ga, angpart)
    # Note: CHECK VARIABLES

    # normalisation of the GTOs
    norma = np.power(np.divide(2, np.pi), 0.75)
    normb = np.zeros(l.size)
    normc = np.zeros(l.size)
    normd = np.zeros(l.size)
    for i in range(l.size):
        normb[i] = np.power(2, l[i] + m[i] + n[i])
        normc[i] = np.power(ga[i],
                            np.divide(np.multiply(2, l[i]) + np.multiply(2, m[i]) + np.multiply(2, n[i]) + 3, 4.0))
        normd[i] = np.multiply(np.multiply(factd(np.multiply(2, l[i]) - 1), factd(np.multiply(2, m[i]) - 1)),
                               factd(np.multiply(2, n[i]) - 1))
    nrm = norma * normb * normc / normd ** 0.5

    # combine the normalisation, MO, and contraction coefficients
    mmod = np.multiply(np.multiply(mos, ci), nrm, dtype=np.double)

    # total number of primitive GTOs
    ncap = l.size
    print('The number of primitive GTOs before reduction: ', str(ncap))

    # generate table with GTO descriptors

    fulltable = np.array(np.vstack((ga, l, m, n, xx, yy, zz)))  # Note: is this table correct?

    ndup = l.size
    # dummy,ipos,irep = unique(full_table,'rows','stable')
    dummy, ipos, irep = uunique(fulltable)
    # ipos = np.sort(ipos)
    # irep = np.sort(irep)
    # print(irep)
    # print(dummy)
    # print(np.size(irep), ncap, ipos.size, "no pasaran, las dimensiones, no pasaran")
    # reduce the matrices to include only the unique cGTOs
    ga = ga[ipos]
    l = np.asarray(l[ipos], dtype=np.int32)  # AMBIGUOUS VARIABLE NAME?
    # print(ipos, irep, "fascistas")
    m = m[ipos]
    n = n[ipos]
    xx = xx[ipos]
    yy = yy[ipos]
    zz = zz[ipos]
    angpart = angpart[ipos]

    # total number of primitive GTOs
    ncap = l.size

    print('The number of primitive GTOs after reduction: ', str(ncap))
    # precomputing the P=0 cases
    nq = q.size
    p0matrix = np.zeros(
        (nq, np.multiply(4, np.max(l)) + 1, np.multiply(4, np.max(l)) + 1, np.multiply(4, np.max(l)) + 1))
    for i in range(np.multiply(4, np.max(l)) + 1):
        for j in range(np.multiply(4, np.max(l)) + 1):
            for k in range(np.multiply(4, np.max(l)) + 1):
                p0matrix[:, i, j, k] = pzerot(i, j, k, q)
    # print("p0 first point:", p0matrix[:, 1, 0, 0])
    m1 = np.asarray(mat[:, 0], dtype=np.int32)
    m2 = np.asarray(mat[:, 1], dtype=np.int32)
    m3 = np.asarray(mat[:, 2], dtype=np.int32)
    m4 = np.asarray(mat[:, 3], dtype=np.int32)

    # Calculation of the prexponential factors Z and Z2 for all N^2 GTO products and N_orb
    z1 = np.zeros((m1.size, ncap, ncap))
    z2 = np.zeros((m1.size, ncap, ncap))

    for i in range(ncap):
        iduplicates = np.argwhere(irep == irep[ipos[i]])
        iduplicates = np.asarray(iduplicates).flatten()
        # print(ipos[i], iduplicates)
        for j in range(ncap):
            jduplicates = np.argwhere(irep == irep[ipos[j]])
            jduplicates = np.asarray(jduplicates).flatten()
            # print(ipos[j], jduplicates)
            for ii in iduplicates:
                for jj in jduplicates:
                    temp1 = total * (mmod[m1, ii] * mmod[m2, jj] + mmod[m1, jj] * mmod[m2, ii])
                    temp2 = mmod[m3, ii] * mmod[m4, jj] + mmod[m3, jj] * mmod[m4, ii]

                    z1[:, i, j] += temp1
                    z2[:, i, j] += temp2

    # Generating the MacMurchie-Davidson coefficients

    dx = mtg(l, xx, ga)  # Note: CHECK AND INSERT FUNCTION "mdtablegen"
    dy = mtg(m, yy, ga)  # Note: s.a.
    dz = mtg(n, zz, ga)  # Note: s.a.

    # Dealing with the orbitals of a similar type together
    dummy, apos, arep = np.unique(angpart, return_index=True, return_inverse=True)  # Note: is "axis=0" correct?

    nnew = apos.size
    print('The number of GTOs after compression is: ', str(nnew))
    # duplicatesang = np.zeros(nnew)
    # for i in range(nnew):
    #     duplicatesang[i] = np.argwhere(arep == arep[apos[i]])

    # total angular momentum
    ll = l + m + n
    ga = ga[apos]
    xx = xx[apos]
    yy = yy[apos]
    zz = zz[apos]
    ll = ll[apos]

    # Hasta aqui funciona

    px = np.zeros((nnew, nnew))
    py = np.zeros((nnew, nnew))
    pz = np.zeros((nnew, nnew))
    e12 = np.zeros((q.size, nnew, nnew))

    # Combined MD coefficients
    ddx = np.zeros((nnew, nnew, 2 * (np.max(l) + 2) + 1, np.max(l) + 1, np.max(l) + 1), dtype=np.float32)
    ddy = np.zeros((nnew, nnew, 2 * (np.max(m) + 2), np.max(m) + 1, np.max(m) + 1), dtype=np.float32)
    ddz = np.zeros((nnew, nnew, 2 * (np.max(n) + 2), np.max(n) + 1, np.max(n) + 1), dtype=np.float32)

    for i in range(nnew):
        iduplicates = np.asarray(np.argwhere(arep == arep[apos[i]]), dtype=np.int32)
        # iduplicates = duplicatesang[i]
        for j in range(nnew):
            jduplicates = np.asarray(np.argwhere(arep == arep[apos[j]]), dtype=np.int32)
            # jduplicates = duplicatesang[j]
            #       various factors
            gap = ga[i] + ga[j]
            px[i, j] = np.divide(np.multiply(ga[i], xx[i]) + np.multiply(ga[j], xx[j]), gap)
            py[i, j] = np.divide(np.multiply(ga[i], yy[i]) + np.multiply(ga[j], yy[j]), gap)
            pz[i, j] = np.divide(np.multiply(ga[i], zz[i]) + np.multiply(ga[j], zz[j]), gap)
            prexp = np.multiply(np.power(np.divide(np.pi, gap), 1.5),
                                np.exp(-np.multiply(np.square(q), np.divide(0.25, gap))))  # Note: shoul "q" be squared?
            #       total pre-exponential factor
            e12[:, i, j] = np.multiply(prexp,
                                       np.exp(-np.multiply(np.divide(np.multiply(ga[i], ga[j]), gap),
                                                           np.square(xx[i] - xx[j]) + np.square(yy[i] - yy[j])
                                                           + np.square(zz[i] - zz[j]))))
            #       combination of the MD coefficients
            iduplicates = iduplicates.flatten()
            jduplicates = jduplicates.flatten()
            for ii in iduplicates:
                for jj in jduplicates:

                    for ls in range(l[ii] + l[jj] + 1):
                        ddx[i, j, ls, l[ii], l[jj]] = dx[ii, jj, ls, l[ii], l[jj]]
                    for ms in range(m[ii] + m[jj] + 1):
                        ddy[i, j, ms, m[ii], m[jj]] = dy[ii, jj, ms, m[ii], m[jj]]
                    for ns in range(n[ii] + n[jj] + 1):
                        ddz[i, j, ns, n[ii], n[jj]] = dz[ii, jj, ns, n[ii], n[jj]]

    # update of the MD coeffs
    dx = ddx
    dy = ddy
    dz = ddz

    # definition of the new size
    ncap = nnew

    # Trying vectorization for the integral_ijkr
    lmax = int(np.max(l))
    lmat = np.zeros((int((lmax + 1) * (lmax + 2) / 2), lmax + 1))
    mmat = np.zeros((int((lmax + 1) * (lmax + 2) / 2), lmax + 1))
    nmat = np.zeros((int((lmax + 1) * (lmax + 2) / 2), lmax + 1))
    maxcount = np.ones(lmax + 1)
    count = 0
    lmaxvec = range(lmax + 1)
    for i in lmaxvec:
        maxcount[i] = (i + 1) * (i + 2) / 2
        count = 0
        for l1 in range(0, i + 1):
            for j in range(0, i - l1 + 1):
                lmat[count, i] = l1
                mmat[count, i] = j
                nmat[count, i] = i - l1 - j
                count += 1

    print('number of values', lmat, count, lmax)

    # Loops over 4xpGTOs
    # Strictly non-diagonal elements with all GTOs being different
    tsi = np.zeros(q.size)
    int_res = np.zeros(q.size)

    for i in range(ncap):
        for j in range(i + 1, ncap):
            for k in range(i + 1, ncap):  # Note: is the range correct? "i+1" or "j+1"? Andres says "i+1".
                for r in range(k + 1, ncap):
                    #               The coordinates of the point H
                    #  print(i,j,k,r)
                    hx = px[k, r] - px[i, j]
                    hy = py[k, r] - py[i, j]
                    hz = pz[k, r] - pz[i, j]
                    h = (hx * hx + hy * hy + hz * hz) ** 0.5
                    # New
                    # print(ll[i], ll[j], ll[k], ll[r])
                    if h < cutoffcentre:

                        #                   compute the zero cases
                        f = intkzero(nq, ll[i], ll[j], ll[k], ll[r], p0matrix, dx, dy, dz, i, j, k, r,
                                     z1, z2, apos, cutoffz, cutoffmd)
                    else:
                        apos2 = np.array(apos, dtype=np.int) + 1
                        #                   computation of the F-integral / sum

                        f = integrals_ijkr.tot_integral_k_ijkr(q, ll[i], ll[j], ll[i], ll[r], hx, hy, hz, h, dx, dy, dz,
                                                               i, j, k, r,
                                                               z1, z2, apos2, cutoffz, cutoffmd)

                    # f = intk_vec(q, angvec1, angvec2, angvec3, angvec4, maxi, maxj, maxk, maxr, ll[i], ll[j], ll[k],
                    #    ll[r], hx, hy, hz, h, dx, dy, dz, i,
                    #  j, k, r, z1, z2, apos, cutoffz, cutoffmd)
                    # f = intk(q, ll[i], ll[j], ll[i], ll[r], hx, hy, hz, h, dx, dy, dz, i, j, k, r,
                    # z1, z2, apos, cutoffz, cutoffmd)

                    #               add to the total intensity
                    tsi += 8 * f * e12[:, i, j] * e12[:, k, r]
    print('hola')
    # diagonal with respect two 1st and 3rd element (k = i)
    for i in range(ncap):
        for j in range(i + 1, ncap):
            for r in range(i + 1, ncap):
                #           coordinates of the point H
                hx = px[i, r] - px[i, j]
                hy = py[i, r] - py[i, j]
                hz = pz[i, r] - pz[i, j]
                h = (hx * hx + hy * hy + hz * hz) ** 0.5
                if h < cutoffcentre:
                    #             computation of the zero cases
                    f = intkzero(nq, ll[i], ll[j], ll[i], ll[r], p0matrix, dx, dy, dz, i, j, i, r,
                                 z1, z2, apos, cutoffz, cutoffmd)
                else:
                    apos2 = np.array(apos, dtype=np.int) + 1
                    #               computation of the F-integral / sum
                    f = integrals_ijkr.tot_integral_k_ijkr(q, ll[i], ll[j], ll[i], ll[r], hx, hy, hz, h, dx, dy, dz, i,
                                                           j, i, r,
                                                           z1, z2, apos2, cutoffz, cutoffmd)
                #           add to the total intensity
                tsi += 4 * f * e12[:, i, j] * e12[:, i, r]

    # diagonal in only one pair (j = i)
    for i in range(ncap):
        for k in range(ncap):
            for r in range(k + 1, ncap):
                #           The coordinates of the point H
                hx = px[k, r] - px[i, i]
                hy = py[k, r] - py[i, i]
                hz = pz[k, r] - pz[i, i]
                h = (hx * hx + hy * hy + hz * hz) ** 0.5
                if h < cutoffcentre:
                    #               compute the zero cases
                    f = intkzero(nq, ll[i], ll[i], ll[k], ll[r], p0matrix, dx, dy, dz, i, i, k, r,
                                 z1, z2, apos, cutoffz, cutoffmd)
                else:
                    apos2 = np.array(apos, dtype=np.int) + 1
                    #               compute the F-integral / sum
                    f = integrals_ijkr.tot_integral_k_ijkr(q, ll[i], ll[i], ll[k], ll[r], hx, hy, hz, h, dx, dy, dz, i,
                                                           i, k, r,
                                                           z1, z2, apos2, cutoffz, cutoffmd)
                #           add to the total intensity
                if i == 1 and k == 2 and r == 9:
                    print(f)
                    print(i, i, k, r, tsi)
                tsi += 4 * f * e12[:, i, i] * e12[:, k, r]

    # diagonal in both pairs without repetition of the pair (j = i, r = k)
    for i in range(ncap):
        for k in range(i + 1, ncap):
            #       The coordinates of the point H
            hx = px[k, k] - px[i, i]
            hy = py[k, k] - py[i, i]
            hz = pz[k, k] - pz[i, i]
            h = (hx * hx + hy * hy + hz * hz) ** 0.5
            if h < cutoffcentre:
                #           compute the zero cases
                f = intkzero(nq, ll[i], ll[i], ll[k], ll[k], p0matrix, dx, dy, dz, i, i, k, k,
                             z1, z2, apos, cutoffz, cutoffmd)
            else:
                apos2 = np.array(apos, dtype=np.int) + 1
                #           compute the F-integral / sum
                f = integrals_ijkr.tot_integral_k_ijkr(q, ll[i], ll[i], ll[k], ll[k], hx, hy, hz, h, dx, dy, dz, i, i,
                                                       k, k,
                                                       z1, z2, apos2, cutoffz, cutoffmd)
            #       add to the total intensity
            tsi += 2 * f * e12[:, i, i] * e12[:, k, k]

    # all GTOs are identical (j = r = k = i)
    for i in range(ncap):
        f = intkzero(nq, ll[i], ll[i], ll[i], ll[i], p0matrix, dx, dy, dz, i, i, i, i,
                     z1, z2, apos, cutoffz, cutoffmd)
        #   add to the total intensity
        tsi += f * e12[:, i, i] * e12[:, i, i]

    print('Maximum intenstiy: ', q[0], np.max(tsi))

# main()
