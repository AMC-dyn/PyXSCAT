import numpy as np


def integral_k_ijkr_0_case(Nq, lmax1, lmax2, lmax3, lmax4, P0matrix, Dx, Dy, Dz,
                           i, j, k, r, Z,Z2, apos, cutOffZ, cutOffMD):
    integ = np.zeros(Nq)

    posI = apos[i]

    for l1 in range(0, lmax1 + 1):
        for m1 in range(0, (lmax1 - l1 + 1)):
            n1 = lmax1 - l1 - m1
            posJ = apos[j]
            for l2 in range(0, lmax2 + 1):
                for m2 in range(0, (lmax2 - l2 + 1)):
                    n2 = lmax2 - l2 - m2
                    Zij = Z[:, posI, posJ]
                    Zij2 = Z2[:, posI, posJ]
                    posK = apos[k]
                    for l3 in range(0, lmax3 + 1):
                        for m3 in range(0, lmax3 - l3 + 1):
                            n3 = lmax3 - l3 - m3
                            posR = apos[r]
                            for l4 in range(0, lmax4 + 1):
                                for m4 in range(0, (lmax4 - l4 + 1)):
                                    n4 = lmax4 - l4 - m4
                                    Zkr = Z[:, posK, posR]
                                    Zkr2 = Z2[:, posK, posR]

                                    Ztot = np.sum(Zij * Zkr2 + Zij2 * Zkr) / 8

                                    if abs(Ztot) < cutOffZ:
                                        posR = posR + 1
                                        continue

                                    for L in range(0, (l1 + l2 + 1)):
                                        MDL = Dx[i, j, L, l1, l2] * Ztot
                                        if MDL == 0:
                                            continue
                                        for M in range(0, m1 + m2 + 1):

                                            MDM = Dy[i, j, M, m1, m2] * MDL
                                            if MDM == 0:
                                                continue

                                            for N in range(0, (n1 + n2 + 1)):
                                                H1 = (-1) ** (L + M + N)
                                                MDN = Dz[i, j, N, n1, n2] * MDM * H1
                                                if MDN == 0:
                                                    continue

                                                for Lp in range(0, (l3 + l4 + 1)):
                                                    MDLp = Dx[k, r, Lp, l3, l4] * MDN
                                                    if MDLp == 0:
                                                        continue
                                                    for Mp in range(0, (m3 + m4 + 1)):
                                                        MDMp = Dy[k, r, Mp, m3, m4] * MDLp
                                                        if MDMp == 0:
                                                            continue
                                                        for Np in range(0, (n3 + n4 + 1)):
                                                            prodD = Dz[k, r, Np, n3, n4] * MDMp

                                                            if abs(prodD) < cutOffMD:
                                                                continue

                                                            f = P0matrix[:, L + Lp, M + Mp, N + Np]
                                                            integ = integ + prodD * f

                                    posR = posR + 1

                            posK = posK + 1

                    posJ = posJ + 1

            posI = posI + 1
    return integ
