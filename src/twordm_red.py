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
    closed = 'none'
    with open('../inputs/abinitio.dat', 'r') as fh2:
        for line in fh2:
            read = False
            running = False
            runlogic = 'none'
            if line.startswith('RUN'):
                running = True
            if running and not line.startswith('RUN'):
                runlogic = line.strip().split()
                running = False

            if not runlogic == 'none' and 'Y' in runlogic:
                if line.startswith('CLOSED'):
                    read = True
                if read and not line.startswith('CLOSED'):
                    closed = line.strip().split()
            else:
                closed = 1
    print('closed orbitals', closed)
    # if 'none' in closed:
    #     # closed = input('Specify the number of closed orbitals \n')
    #     closed = 0
    countcivs=0
    closed = int(closed)
    with open('molpro.pun', 'r') as fh:
        for line in fh:
            if not line[0].isupper() and not line[0] == '*' and not line[0] == '?' and not line[0] == '-':
                NN = line.strip().split()
                if "." not in NN[0]:
                    civs.append(NN[1:])
                    confs.append(
                        'ab' * closed + NN[0].replace("0", "00").replace("a", "a0").replace("b", "0b").replace("2",
                                                                                                        "ab"))

                    countcivs += 1
                else:

                    civs[countcivs-1].extend(NN)
                if len(NN) != 1:
                    count += 1
    # up to here in python

    # twordm = np.zeros((int(len(confs[0]) / 2), int(len(confs[0]) / 2), int(len(confs[0]) / 2), int(len(confs[0]) / 2)))
    # nc1 = 0
    # ep = 1
    print('read complete with ', np.size(confs) * np.size(confs))
    # ep2 = np.ones((len(confs), len(confs)))
    # for c1 in confs:
    #     nc2 = 0
    #     for c2 in confs:
    #
    #         mat1 = [i for i in range(len(c1)) if c1[i] != "0"]
    #         mat2 = [i for i in range(len(c1)) if c2[i] != "0"]
    #
    #         for i in range(len(mat1)):
    #
    #             if mat1[i] != mat2[i]:
    #                 for j in range(len(mat2)):
    #                     if mat1[i] == mat2[j]:
    #                         ep2[nc1, nc2] = -ep2[nc1, nc2]
    #
    #         nc2 += 1
    #     nc1 += 1
    # #
    # ndiff = np.ones((np.size(confs), np.size(confs)))
    # nc1 = 0
    # print('maybe this part takes forever')
    #
    # for c1 in confs:
    #     nc2 = 0
    #     for c2 in confs:
    #         ndiff[nc1, nc2] = sum(1 for a, b in zip(c1, c2) if a != b)
    #
    #         nc2 = nc2 + 1
    #     nc1 = nc1 + 1
    # print('number of differences', np.shape(ndiff))
    # print(ep2)
   # print(civs[0])
    return civs, confs
