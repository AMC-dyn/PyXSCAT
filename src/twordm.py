# -*- coding: utf-8 -*-
"""
Created on Sun Jul 11 17:06:33 2021

@author: AndresMoreno
"""

import numpy as np


def twordmconst():
    count = 0
    confs = []
    civs = []

    with open('molpro.pun', 'r') as fh:
        for line in fh:
            if not line[0].isupper() and not line[0] == '*' and not line[0] == '?' and not line[0] == '-':
                NN = line.strip().split()

                civs.append(NN[1:])
                confs.append(NN[0].replace("0", "00").replace("a", "a0").replace("b", "0b").replace("2", "ab"))
                if len(NN) != 1:
                    count += 1

    twordm = np.zeros((int(len(confs[0]) / 2), int(len(confs[0]) / 2), int(len(confs[0]) / 2), int(len(confs[0]) / 2)))
    nc1 = 0
    ep = 1
    ep2 = np.ones((len(confs), len(confs)))
    for c1 in confs:
        nc2 = 0
        for c2 in confs:

            mat1 = [i for i in range(len(c1)) if c1[i] != "0"]
            mat2 = [i for i in range(len(c1)) if c2[i] != "0"]

            for i in range(len(mat1)):

                if mat1[i] != mat2[i]:
                    for j in range(len(mat2)):
                        if mat1[i] == mat2[j]:
                            ep2[nc1, nc2] = -ep2[nc1, nc2]

            nc2 += 1
        nc1 += 1

    nc1 = 0
    for c1 in confs:
        nc2 = 0
        for c2 in confs:

            ndiff = sum(1 for a, b in zip(c1, c2) if a != b)
            ep = ep2[nc1, nc2]
            if ndiff != 0:

                sdef = False
                rdef = False
                pdef = False
                qdef = False

                diffs1 = []
                diffs2 = []
                spin1 = []
                spin2 = []
                for n in range(0, len(c1)):
                    if c1[n] != c2[n]:
                        if c1[n] != "0":
                            diffs1.append(int(np.round((n + 1) / 2 + 0.1)))
                            spin1.append(c1[n])
                        elif c2[n] != "0":
                            diffs2.append(int(np.round((n + 1) / 2 + 0.1)))
                            spin2.append(c2[n])

                if ndiff == 4:

                    sorb = diffs2[0]
                    qorb = diffs2[1]
                    porb = diffs1[1]
                    rorb = diffs1[0]

                    eg = 1.00

                    if spin2[0] == spin1[0] and spin2[1] == spin1[1]:
                        twordm[porb - 1, rorb - 1, sorb - 1, qorb - 1] += float(civs[nc1][0]) * float(
                            civs[nc2][0]) * ep * eg

                    sorb = diffs2[0]
                    qorb = diffs2[1]
                    rorb = diffs1[1]
                    porb = diffs1[0]
                    eg = -1.00

                    if spin2[0] == spin1[1] and spin2[1] == spin1[0]:
                        twordm[porb - 1, rorb - 1, sorb - 1, qorb - 1] += float(civs[nc1][0]) * float(
                            civs[nc2][0]) * ep * eg

                    qorb = diffs2[0]
                    sorb = diffs2[1]
                    porb = diffs1[1]
                    rorb = diffs1[0]
                    eg = -1.00

                    if spin1[0] == spin2[1] and spin2[0] == spin1[1]:
                        twordm[porb - 1, rorb - 1, sorb - 1, qorb - 1] += float(civs[nc1][0]) * float(
                            civs[nc2][0]) * ep * eg

                    qorb = diffs2[0]
                    sorb = diffs2[1]
                    rorb = diffs1[1]
                    porb = diffs1[0]
                    eg = 1.00

                    if spin1[1] == spin2[1] and spin2[0] == spin1[0]:
                        twordm[porb - 1, rorb - 1, sorb - 1, qorb - 1] += float(civs[nc1][0]) * float(
                            civs[nc2][0]) * ep * eg


                elif ndiff == 2:

                    qorb = diffs2[0]
                    porb = diffs1[0]
                    eg = 1.00

                    for i in range(1, int(len(c2) + 1)):

                        if c2[i - 1] != 0:
                            sdef = True

                            sorb = int(round((i) / 2 + 0.1))
                            # print("proof:", sorb,porb)
                            spins = c2[i - 1]

                        if c1[i - 1] != "0":
                            # print("nayy")
                            rdef = True
                            rorb = int(round((i) / 2 + 0.1))
                            spinr = c1[i - 1]

                        if sdef and rdef and spins == spinr and spin2[0] == spin1[0]:
                            sdef = False
                            rdef = False
                            pdef = False
                            qdef = False

                            twordm[porb - 1, rorb - 1, sorb - 1, qorb - 1] += float(civs[nc1][0]) * float(
                                civs[nc2][0]) * ep * eg



                        else:
                            sdef = False
                            rdef = False
                            pdef = False
                            qdef = False

                    qorb = diffs2[0]
                    rorb = diffs1[0]
                    eg = -1.00

                    for i in range(1, len(c2) + 1):

                        if c2[i - 1] != "0":
                            sdef = True
                            sorb = round((i) / 2 + 0.1)
                            spins = c2[i - 1]

                        if c1[i - 1] != "0":
                            pdef = True
                            porb = round((i) / 2 + 0.1)
                            spinp = c1[i - 1]

                        if sdef and pdef and spin1[0] == spins and spin2[0] == spinp:
                            sdef = False
                            rdef = False
                            pdef = False
                            qdef = False
                            twordm[porb - 1, rorb - 1, sorb - 1, qorb - 1] += float(civs[nc1][0]) * float(
                                civs[nc2][0]) * ep * eg


                        else:
                            sdef = False
                            rdef = False
                            pdef = False
                            qdef = False

                    sorb = diffs2[0]
                    porb = diffs1[0]
                    eg = -1.00

                    for i in range(1, len(c2) + 1):

                        if c2[i - 1] != "0":
                            qdef = True
                            qorb = round((i) / 2 + 0.1)
                            spinq = c2[i - 1]

                        if c1[i - 1] != "0":
                            rdef = True
                            rorb = round((i) / 2 + 0.1)
                            spinr = c1[i - 1]

                        if rdef and qdef and spins == spinr and spin1[0] == spin2[0]:
                            sdef = False
                            rdef = False
                            pdef = False
                            qdef = False
                            twordm[porb - 1, rorb - 1, sorb - 1, qorb - 1] += float(civs[nc1][0]) * float(
                                civs[nc2][0]) * ep * eg


                        else:
                            sdef = False
                            rdef = False
                            pdef = False
                            qdef = False

                    sorb = diffs2[0]
                    rorb = diffs1[0]
                    eg = 1.00

                    for i in range(1, len(c2) + 1):
                        if c2[i - 1] != "0":
                            qdef = True
                            qorb = round((i) / 2 + 0.1)
                            spinq = c2[i - 1]

                        if c1[i - 1] != "0":
                            pdef = True
                            porb = round((i) / 2 + 0.1)
                            spinp = c1[i - 1]

                        if qdef and pdef and spinq == spinp and spin2[0] == spin1[0]:
                            qdef = False
                            pdef = False

                            twordm[porb - 1, rorb - 1, sorb - 1, qorb - 1] += float(civs[nc1][0]) * float(
                                civs[nc2][0]) * ep * eg

                        else:

                            pdef = False
                            qdef = False
            elif ndiff == 0:

                ep = 1

                for i1 in range(1, len(c1) + 1):

                    for i2 in range(1, len(c2) + 1):

                        if i1 != i2:

                            if c1[i1 - 1] != "0" and c1[i2 - 1] != "0":

                                sorb = int(np.round(i1 / 2 + 0.1))
                                qorb = int(np.round(i2 / 2 + 0.1))

                                if c1[i1 - 1] == c2[i2 - 1]:
                                    porb = sorb
                                    rorb = qorb

                                    eg = -1.00

                                    twordm[porb - 1, rorb - 1, sorb - 1, qorb - 1] += float(civs[nc1][0]) * float(
                                        civs[nc2][0]) * ep * eg

                                porb = qorb
                                rorb = sorb
                                eg = 1.00

                                twordm[porb - 1, rorb - 1, sorb - 1, qorb - 1] += float(civs[nc1][0]) * float(
                                    civs[nc2][0]) * ep * eg

            nc2 += 1
        nc1 += 1
    cutoff = 1e-9
    mat = []
    total = []
    for p in range(int(len(confs[0]) / 2)):
        for q in range(int(len(confs[0]) / 2)):
            for r in range(int(len(confs[0]) / 2)):
                for s in range(int(len(confs[0]) / 2)):
                    if abs(twordm[p, q, r, s]) >= cutoff:
                        mat.append([p, s, q, r])
                        total.append(twordm[p, q, r, s])
    mat = np.asarray(mat)
    total = np.asarray(total)
    f = open('twordm.dat', 'w')
    for i in range(np.size(total)):
        f.write("{0} {1}\n".format(mat[i, :].flatten(), str(total[i])))
    f.close()
    mat = np.asarray(mat, dtype=np.int64)
    np.savetxt('barrileros.dat', np.vstack((np.transpose(mat), total)))
    # np.savetxt('patron.dat', total)

    return np.asarray(mat, dtype=np.int64), np.asarray(total, dtype=np.float64)
