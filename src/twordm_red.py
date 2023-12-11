# -*- coding: utf-8 -*-
"""
Created on 19 July 2023

@authors: Mats Simmermacher and Andres Moreno Carrascosa
"""

import numpy as np
import sys
import symmetry as sym

def twordmconst(closed,fileread,molpro,molpro2022):

    # check of the closed orbitals:

    if isinstance(closed, list):
        ncirr = len(closed)
        nclosed = sum(closed)
        if ncirr == 1:
            closed = closed[0]
    else:
        ncirr = 1
        nclosed = closed

    # extract the number of states and the number of irreducible representations:
    states1 = []
    states2 = []
    if not molpro2022:

        with open(fileread, 'r') as fh:
            for line in fh:
                if line[0].isupper():
                    NN = line.strip().split()
                    if NN[0] == 'MCSCF' and ( NN[1] == 'STATE' or NN[1] == 'STATES' ):
                        states1.append(float(NN[2]))
                elif not line[0] == '*' and not line[0] == '?' and not line[0] == '-':
                    firstline = line.strip().split()
                    break
        states1 = np.unique(states1)
        nstates1 = len(states1)
        for state in states1:
            states2.append(round(10 * (state - np.floor(state))))
        states2 = np.unique(states2)
        nstates2 = len(states2)
        nirr = -1
        for element in firstline:
            nirr = nirr + 1
            if '.' in element:
                break

        # checks whether the specification of the closed orbitals matches the number of irreps:

        if ncirr != nirr:
            sys.exit(print('Error: mismatch between the closed orbital entry and the number of irreps in the CI vector.', nirr, ' entries are required.'))
        print('total number of closed orbitals: ', nclosed)
        print('number of irreducible representations in CI vector: ', nirr)
        print('total number of states: ', nstates1)
        print('number of irreducible representations of states: ', nstates2)

        # read the CI vector

        # if no symmetry is used or present:

        if nirr == 1:
            count = 0
            countcivs = 0
            civs = []
            confs = []
            with open(fileread, 'r') as fh:
                for line in fh:
                    if not line[0].isupper() and not line[0] == '*' and not line[0] == '?' and not line[0] == '-':
                        NN = line.strip().split()
                        if '.' not in NN[0]:
                            sdet = ['c' * closed]
                            sdet.append(NN[0])
                            sdet = ''.join(sdet)
                            sdet = sdet.replace('0', '00').replace('a', 'a0').replace('b', '0b').replace('2', 'ab')
                            sdet = sdet.replace('c', 'ab')
                            confs.append(sdet)
                            civs.append(NN[1:])
                            countcivs += 1
                        else:
                            civs[countcivs-1].extend(NN)
                        if len(NN) != 1:
                            count += 1
        # if states are only present in one irreducible representation:
        elif nstates2 == 1:
            count = 0
            countcivs = 0
            civs = []
            confs = []
            with open(fileread, 'r') as fh:
                for line in fh:
                    if not line[0].isupper() and not line[0] == '*' and not line[0] == '?' and not line[0] == '-':
                        NN = line.strip().split()
                        if '.' not in NN[0]:
                            sdet = []
                            for irr in range(nirr):
                                sdet.append('c' * closed[irr])
                                sdet.append(NN[irr])
                            sdet = ''.join(sdet)
                            sdet = sdet.replace('0', '00').replace('a', 'a0').replace('b', '0b').replace('2', 'ab')
                            sdet = sdet.replace('c', 'ab')
                            confs.append(sdet)
                            civs.append(NN[nirr:])
                            countcivs += 1
                        else:
                            civs[countcivs-1].extend(NN)
                        if len(NN) != 1:
                            count += 1
        # if states are present in more than one irreducible representation:
        else:
            count = 0
            countcivs = 0
            civs = np.empty((nstates2,0)).tolist()
            confs = np.empty((nstates2,0)).tolist()
            # irrstr = np.empty((nstates2,0)).tolist()
            with open(fileread, 'r') as fh:
                for line in fh:
                    if not line[0].isupper() and not line[0] == '*' and not line[0] == '?' and not line[0] == '-':
                        NN = line.strip().split()
                        if '.' not in NN[0]:
                            irrep = sym.findirrep(nirr,NN[:nirr])
                            # irrstr[irrep-1].append(irrep)
                            sdet = []
                            for irr in range(nirr):
                                sdet.append('c' * closed[irr])
                                sdet.append(NN[irr])
                            sdet = ''.join(sdet)
                            sdet = sdet.replace('0', '00').replace('a', 'a0').replace('b', '0b').replace('2', 'ab')
                            sdet = sdet.replace('c', 'ab')
                            confs[irrep-1].append(sdet)
                            civs[irrep-1].append(NN[nirr:])
                            countcivs += 1
                        else:
                            civs[irrep-1,countcivs-1].extend(NN)
                        if len(NN) != 1:
                            count += 1

    elif molpro2022:
        count = 0
        countcivs = 0
        civs = []
        confs = []
        with open(fileread, 'r') as fh:
            for line in fh:
                if "CI Coefficients of symmetry 1" in line:
                    break
            fh.readline()
            fh.readline()
            for line in fh:
                if "Energy:" in line or line.isspace():
                    break

                NN = line.strip().split()
                sdet = ['c' * closed]
                sdet.append(NN[0])
                sdet = ''.join(sdet)
                sdet = sdet.replace('0', '00').replace('a', 'a0').replace('b', '0b').replace('2', 'ab')
                sdet = sdet.replace('c', 'ab')
                confs.append(sdet)
                civs.append(NN[1:])
                countcivs += 1





    print('read complete with ', np.size(confs) * np.size(confs))
    
    return civs, confs
