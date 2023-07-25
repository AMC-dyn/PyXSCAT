import itertools as it
import re
import sys

import h5py
import numpy as np

# converter


def csf_to_slater_basis_conversion(csfs, civs, S):
    sds = []
    dictionary = {}
    civecinsd = []
    thresh = 1e-8
    # convert csfs into sds, and then multiply the coeff in csf basis by conversion
    length = len(csfs)
    for i in range(len(csfs)):
        if i % 1000 == 0:
            print(f"{i}/{length}")
        if abs(civs[i]) < thresh:
            continue

        conv = csf2sd(csfs[i], S)

        for j in conv:
            #  if abs(j[0]) < thresh:
                #  continue
            # if this is slow, I know a way to fix it
            try:
                dictionary[j[1]] += j[0] * float(civs[i])
            except KeyError:
                dictionary[j[1]] = j[0] * float(civs[i])

            #  if j[1] not in sds:
            #  sds.append(j[1])
            #  civecinsd.append(0.)
            #  indx = sds.index(j[1])
            #  civecinsd[indx] += j[0] * float(civs[i])


#
#  for i in range(len(civecinsd)):
#  civecinsd[i] = str(civecinsd[i])

    sds = list(dictionary.keys())
    print(len(sds))
    civecinsd = [dictionary[key] for key in sds]
    print(np.sum(np.array(civecinsd)**2))

    return sds, civecinsd


def create_sym_csf(string, frozen, inact, ras2, second):
    no_irreps = len(frozen)
    csf = ''
    j = 0
    for i in range(no_irreps):
        #  if frozen[i]+inact[i]+ras2[i]+second[i] < 1:
        #  continue
        csf += frozen[i] * '2'
        csf += inact[i] * '2'
        if ras2[i] > 0:
            csf += string[j]
            j += 1
        #  csf += second[i]*'0'
    return csf


def csf2sd(csfvec, m_s):
    ''' takes in a string of form (e.g.) "22ud00", and returns a list containing all of the slater determinents used in expanding that CSF.

    The list is in the form 
    det[i][0] = value of expansion coefficient
    det[i][1] = string of expansion, using standard 2,a,b,0 notation
    '''

    det = []
    no_mo = len(csfvec)

    no_d = csfvec.count('d')
    no_u = csfvec.count('u')

    # formulate  Paldus vectors
    a = np.zeros(no_mo)
    b = np.zeros(no_mo)
    c = np.zeros(no_mo)

    curr_orb = csfvec[0]

    if curr_orb == '0':
        c[0] += 1
    elif curr_orb == 'u':
        b[0] += 1
    elif curr_orb == 'd':
        a[0] += 1
        b[0] -= 1
        c[0] += 1
    elif curr_orb == '2':
        a[0] += 1
    else:
        print(csfvec)
        print('ERROR - not a CSF')
        exit

    for i in range(1, no_mo):
        a[i] = a[i - 1]
        b[i] = b[i - 1]
        c[i] = c[i - 1]
        curr_orb = csfvec[i]
        if curr_orb == '0':
            c[i] += 1
        elif curr_orb == 'u':
            b[i] += 1
        elif curr_orb == 'd':
            a[i] += 1
            b[i] -= 1
            c[i] += 1
        elif curr_orb == '2':
            a[i] += 1
        else:
            print('ERROR - not a CSF')
            exit

    no_a = int((no_u + no_d + 2 * m_s) / 2)

    slat_det = [None] * no_mo

    for alpha_somos in it.combinations(range(0, no_u + no_d), no_a):
        coeff = 1
        i_s = 0
        i_a = 0
        i_b = 0
        for i in range(no_mo):
            curr_orb = csfvec[i]
            if curr_orb == '0':
                slat_det[i] = '0'
            elif curr_orb == 'u':
                if i_s in alpha_somos:
                    slat_det[i] = 'a'
                    coeff *= (a[i] + b[i] - i_b)
                    i_a += 1
                else:
                    slat_det[i] = 'b'
                    coeff *= (a[i] + b[i] - i_a)
                    i_b += 1
                coeff /= b[i]
                i_s += 1
            elif curr_orb == 'd':
                if i_s in alpha_somos:
                    slat_det[i] = 'a'
                    coeff *= (i_b - a[i] + 1)
                    i_a += 1
                    if b[i] % 2 == 0:
                        coeff *= -1
                else:
                    slat_det[i] = 'b'
                    coeff *= (i_a - a[i] + 1)
                    i_b += 1
                    if b[i] % 2 != 0:
                        coeff *= -1
                coeff /= (b[i] + 2)
                i_s += 1
            elif curr_orb == '2':
                slat_det[i] = '2'
                if b[i] % 2 != 0:
                    coeff *= -1
                i_a += 1
                i_b += 1

        if coeff != 0:

            det.append(
                [np.sign(coeff) * np.sqrt(np.abs(coeff)), ''.join(slat_det)])

    return det


def civec2bin(a):
    binary = ''
    for i in a:
        if i == '2':
            binary += '11'
        elif i == 'u':
            binary += '10'
        elif i == 'd':
            binary += '01'
        elif i == '0':
            binary += '00'
        else:
            print("CSF doesn't look correct - check the output")
    return int(binary, 2)


def get_civecs_in_csfs(filename, caspt2):
    with open(filename) as f:
        civs = []
        for line in f:
            if re.search('.*Symmetry species.*', line):
                no_irreps = len(line.split()) - 2
                if no_irreps > 1:
                    print('Using symmetry in the calculations')
                    sym = True
                    line = f.readline()
                    irrep_numbers = [int(a) for a in line.split()[2:]]
                    no_irreps = len(irrep_numbers)
                else:
                    sym = False
                break
        for line in f:
            if re.search('.*Number of closed shell electrons.*', line):
                no_closed = int(line.split()[5]) // 2
                print('no_closed', no_closed)
                break
        for line in f:
            if re.search('.*Spin quantum number.*', line):
                S = int(float(line.split()[3]))
                break

        if sym:
            for line in f:
                if re.search('.*Frozen orbitals.*', line, re.IGNORECASE):
                    no_frozen = [int(x) for x in line.split()[2:]]
                    print('Frozen', no_frozen)
                    break
            for line in f:
                if re.search('.*Inactive orbitals.*', line, re.IGNORECASE):
                    no_inact = [int(x) for x in line.split()[2:]]
                    print('Inactive', no_inact)
                    break
            for line in f:
                if re.search('.*RAS2 orbitals.*', line, re.IGNORECASE):
                    no_ras2 = [int(x) for x in line.split()[2:]]
                    print('Ras2', no_ras2)
                    no_act_irreps = sum(i > 0 for i in no_ras2)
                    print(no_act_irreps, ' active irreps')
                    break
            for line in f:
                if re.search('.*Secondary orbitals.*', line, re.IGNORECASE):
                    no_second = [int(x) for x in line.split()[2:]]
                    print('Secondary', no_second)
                    break
            # Check to see if they add up
            check_var = [
                sum(i) for i in zip(no_frozen, no_inact, no_ras2, no_second)
            ]
            print(check_var)
            if check_var != irrep_numbers:
                print('Help! Things dont add up')
                print(check_var)
                print(irrep_numbers)
                sys.exit()
            for line in f:
                if re.search('.*Number of CSFs.*', line, re.IGNORECASE):
                    no_csfs = int(line.split()[3])
                    break
            if caspt2:  # read a bit further on in the file to ensure CASPT2
                for line in f:
                    if re.search('.*Number of CI roots used.*', line,
                                 re.IGNORECASE):
                        no_roots = int(line.split()[5])
                        break
                for root in range(no_roots):
                    count = 0
                    while True:
                        line = f.readline()
                        if not line:
                            break
                        if re.search(
                                '.*The CI coefficients for the MIXED state nr.',
                                line):
                            csfs = []
                            data = []
                            f.readline()  # line of '-----...'
                            f.readline(
                            )  # line of 'CI COEFFICIENTS LARGER THAN...'
                            f.readline()  # line of 'Occupation..'
                            f.readline()  # line of 'of open shells'
                            f.readline()  # line of 'SGUGA info...'
                            f.readline()  # line of 'Conf  SGUGA...'
                            # beter would be to use the index, but for small active space cases this is probably fine.
                            for j in range(no_csfs):
                                intvar = f.readline()[30:].split()
                                csfs.append(
                                    create_sym_csf(intvar[:-2], no_frozen,
                                                   no_inact, no_ras2,
                                                   no_second))
                                data.append(intvar[-2])
                                count += float(intvar[-2])**2
                            civs.append(data)
                    print('Normalisation for root ', root, ' = ', count)

            else:
                for line in f:
                    if re.search('.*Number of root\(s\) required.*', line,
                                 re.IGNORECASE):
                        no_roots = int(line.split()[4])
                        break
                for root in range(no_roots):
                    count = 0
                    while True:
                        line = f.readline()
                        if not line:
                            break
                        if re.search(
                                '.*printout of CI-coefficients larger than.*',
                                line):
                            csfs = []
                            data = []
                            energy = f.readline().split()[1]
                            f.readline()
                            for j in range(no_csfs):
                                intvar = f.readline().split()[1:2 +
                                                              no_act_irreps]
                                csfs.append(
                                    create_sym_csf(intvar[:-1], no_frozen,
                                                   no_inact, no_ras2,
                                                   no_second))
                                data.append(intvar[-1])
                                count += float(intvar[-1])**2
                            civs.append(data)
                    print('Normalisation for root ', root, ' = ', count)
        else:
            for line in f:
                if re.search('.*Number of CSFs.*', line, re.IGNORECASE):
                    no_csfs = int(line.split()[3])
                    break
            if caspt2:  # read a bit further on in the file to ensure CASPT2
                print('Searching for number of roots in caspt2 calculation')
                for line in f:
                    no_roots = 1
                    if re.search('.*Number of CI roots used.*', line,
                                 re.IGNORECASE):
                        no_roots = int(line.split()[5])
                        break
                if no_roots == 1:
                    f.seek(0)

                for root in range(no_roots):
                    count = 0
                    while True:
                        line = f.readline()
                        if not line:
                            break
                        if re.search(
                                '.*The CI coefficients for the MIXED state nr.',
                                line):
                            print(f'Found state {root}')
                            csfs = []
                            data = []
                            f.readline()  # line of '-----...'
                            f.readline(
                            )  # line of 'CI COEFFICIENTS LARGER THAN...'
                            f.readline()  # line of 'Occupation..'
                            f.readline()  # line of 'of open shells'
                            f.readline()  # line of 'SGUGA info...'
                            f.readline()  # line of 'Conf  SGUGA...'
                            # beter would be to use the index, but for small active space cases this is probably fine.
                            for j in range(no_csfs):
                                intvar = f.readline()[30:].split()
                                csfs.append(intvar[0])
                                data.append(intvar[1])
                                count += float(intvar[1])**2
                            civs.append(data)
                    print('Normalisation for root ', root, ' = ', count)

            else:
                for line in f:
                    if re.search('.*Number of root\(s\) required.*', line,
                                 re.IGNORECASE):
                        no_roots = int(line.split()[4])
                        break
                for root in range(no_roots):
                    count = 0
                    while True:
                        line = f.readline()
                        if not line:
                            break
                        if re.search(
                                '.*printout of CI-coefficients larger than',
                                line):
                            csfs = []
                            data = []
                            energy = f.readline().split()[1]
                            f.readline()
                            for j in range(no_csfs):
                                intvar = f.readline().split()[1:3]
                                csfs.append(intvar[0])
                                data.append(intvar[1])
                                count += float(intvar[1])**2
                            civs.append(data)
                    print('Normalisation for root ', root, ' = ', count)

    # Very important! leave
    if sym:
        no_closed = 0

    return no_roots, S, no_closed, csfs, civs


def sd_formatter(sds, no_closed):
    for i in range(len(sds)):
        sds[i] = no_closed * 'ab' + sds[i].replace("0", "00").replace(
            "a", "a0").replace("b", "0b").replace("2", "ab")
    return sds


def transposer(array):

    return list(zip(*array))


def get_civs_and_confs(logfile, caspt2, extra_precision):
    no_roots, S, no_closed, csfs, civs = get_civecs_in_csfs(logfile, caspt2)
    if extra_precision:
        print('Using extra precision from h5 files')
        if caspt2:
            f = h5py.File('molcas.caspt2.h5', 'r')
            civs = list(f['CI_VECTORS'])
        else:
            f = h5py.File('molcas.rasscf.h5', 'r')
            #  print(civs)
            civs = list(f['CI_VECTORS'])
            #  print(civs)

    sds = []
    newcivec = []
    for i in range(no_roots):
        S = 0.5
        a, b = csf_to_slater_basis_conversion(csfs, civs[i], S)
        sds = a
        newcivec.append(b)

    civs = transposer(newcivec)
    confs = sd_formatter(sds, no_closed)
    return civs, confs


def main():
    filename = 'molcas.log'
    caspt2 = True

    no_roots, S, no_closed, csfs, civs = get_civecs_in_csfs(filename, caspt2)

    sds = []
    newcivec = []
    for i in range(no_roots):
        a, b = csf_to_slater_basis_conversion(csfs, civs[i], S)
        sds = a
        newcivec.append(b)

    sds = sd_formatter(sds, no_closed)

    print(sds, newcivec)


if __name__ == "__main__":
    main()
