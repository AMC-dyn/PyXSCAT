import numpy as np


def BesselDeriv(LL, MM, NN, a, b, c, LLmax):
    TempBesselPre = np.zeros((LLmax + 1))
    for ii in range(0, LL + 1):
        C1 = a[LL, ii]
        if C1 == 0:
            continue
        for jj in range(0, MM + 1):
            Ct2 = b[MM, jj]
            if Ct2 == 0:
                continue
            C2 = C1 * Ct2
            for kk in range(0, NN + 1):
                Ct3 = c[NN, kk]
                if Ct3 == 0:
                    continue
                C3 = C2 * Ct3
                hOrder = np.ceil((LL + MM + NN - ii - jj - kk) / 2) + ii + jj + kk
                hOrder = np.asarray(hOrder, dtype=np.int32)

                TempBesselPre[hOrder] += C3

    return TempBesselPre


def integral_k_ijkr(mu, lmax1, lmax2, lmax3, lmax4, Hx, Hy, Hz, H,
                    Dx, Dy, Dz, i, j, k, r, Z, Z2, apos, cutOffZ, cutOffMD):
    LLmax = lmax1 + lmax2 + lmax3 + lmax4

    Z = np.asarray(Z)
    Z2 = np.asarray(Z2)
    dxij = Dx[i, j, :, :, :]
    dxkr = Dx[k, r, :, :, :]
    dyij = Dy[i, j, :, :, :]
    dykr = Dy[k, r, :, :, :]
    dzij = Dz[i, j, :, :, :]
    dzkr = Dz[k, r, :, :, :]

    a = np.zeros((LLmax + 1, LLmax + 1))
    a[0, 0] = 1

    if LLmax > 0:
        a[1, 1] = -Hx

    for L in range(1, LLmax):
        a[L + 1, 0] = -(L) * a[L - 1, 0]

        for ka in range(0, LLmax):
            a[L + 1, ka + 1] = -Hx * a[L, ka] - L * a[L - 1, ka + 1]

    b = np.zeros((LLmax + 1, LLmax + 1))
    b[0, 0] = 1
    if LLmax > 0:
        b[1, 1] = -Hy

    for M in range(1, LLmax):
        b[M + 1, 0] = -(M) * b[M - 1, 0]

        for ka in range(0, LLmax):
            b[M + 1, ka + 1] = -Hy * b[M, ka] - (M) * b[M - 1, ka + 1]

    c = np.zeros((LLmax + 1, LLmax + 1))
    c[0, 0] = 1
    if LLmax > 0:
        c[1, 1] = -Hz

    for N in range(1, LLmax):
        c[N + 1, 0] = -(N) * c[N - 1, 0]

        for ka in range(0, LLmax):
            c[N + 1, ka + 1] = -Hz * c[N, ka] - (N) * c[N - 1, ka + 1]

    h_saved = np.zeros((LLmax + 1))

    h_pre2 = np.zeros((LLmax + 1, LLmax + 1, LLmax + 1, LLmax + 1))

    for l in range(0, LLmax + 1):
        for m in range(0, LLmax + 1 - l):
            for n in range(0, LLmax + 1 - l - m):
                Temp_h_pre = BesselDeriv(l, m, n, a, b, c, LLmax)
                h_pre2[:, l, m, n] = Temp_h_pre
    posI = apos[i]
    count = 0
    for l1 in range(0, lmax1 + 1):
        for m1 in range(0, (lmax1 + 1 - l1)):
            n1 = lmax1 - l1 - m1
            posJ = apos[j]

            for l2 in range(0, lmax2 + 1):
                for m2 in range(0, (lmax2 + 1 - l2)):
                    n2 = lmax2 - l2 - m2
                    Zij = Z[:, posI, posJ]
                    Zij2 = Z2[:, posI, posJ]

                    posK = apos[k]

                    for l3 in range(0, lmax3 + 1):
                        for m3 in range(0, (lmax3 + 1 - l3)):
                            n3 = lmax3 - l3 - m3

                            posR = apos[r]

                            for l4 in range(0, lmax4 + 1):
                                for m4 in range(0, (lmax4 + 1 - l4)):

                                    n4 = lmax4 - l4 - m4
                                    Zkr = Z[:, posK, posR]
                                    Zkr2 = Z2[:, posK, posR]
                                    varztot = (Zij * Zkr2 + Zij2 * Zkr) * 0.125

                                    Ztot = np.sum(varztot)
                                    count += 1
                                    if i == 0 and j == 1 and k == 1 and r == 16:
                                        print(posI,posJ,posK,posR,Ztot)
                                        print(count)

                                    if abs(Ztot) < cutOffZ:
                                        posR = posR + 1
                                        continue
                                    for L in range(0, (l1 + l2) + 1):
                                        # MDL = Dx[i, j, L, l1, l2] * Ztot
                                        MDL = dxij[L, l1, l2] * Ztot

                                        if MDL == 0:
                                            continue

                                        for M in range(0, (m1 + m2 + 1)):
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
                                                    LL = L + Lp
                                                    for Mp in range(0, (m3 + m4 + 1)):
                                                        MDMp = Dy[k, r, Mp, m3, m4] * MDLp

                                                        if MDMp == 0:
                                                            continue
                                                        MM = M + Mp
                                                        for Np in range(0, (n3 + n4 + 1)):

                                                            prodD = Dz[k, r, Np, n3, n4] * MDMp

                                                            if abs(prodD) < cutOffMD:
                                                                continue

                                                            NN = N + Np
                                                            h_saved += h_pre2[:, LL, MM, NN] * prodD
                                    posR = posR + 1
                            posK = posK + 1
                    posJ = posJ + 1
            posI = posI + 1

    Pmu = H * mu

    h_0 = np.sin(Pmu) / Pmu
    coeff = h_saved[0]
    h_sum = h_0 * coeff

    if LLmax == 1:
        h_1 = (np.sin(Pmu) / Pmu ** 2 - np.cos(Pmu) / Pmu) * mu / H
        coeff = h_saved[0]
        h_sum = h_sum + coeff * h_1
    elif LLmax > 1:
        muOH = mu / H
        h_1 = (np.sin(Pmu) / Pmu ** 2 - np.cos(Pmu) / Pmu) * muOH
        coeff = h_saved[1]
        h_sum = h_sum + coeff * h_1
        for ra in range(2, LLmax):
            coeff = h_saved[ra]
            h_r = ((2 * ra - 1) / Pmu * h_1 - h_0 * muOH) * muOH
            h_sum = h_sum + h_r * coeff
            h_0 = h_1
            h_1 = h_r
    if i == 0 and j == 1 and k == 1 and r == 16:
        print('comparison with integral_ijkr',h_sum)
        print(count)
    return h_sum
