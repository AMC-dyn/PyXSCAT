import itertools as it
import numpy as np



def csf2sd(csfvec, m_s):
    ''' takes in a string of form (e.g.) "22ud00", and returns a list containing all of the slater determinents used in expanding that CSF.

    The list is in the form
    det[i][0] = sign of expansion coefficient
    det[i][1] = value of expansion coefficient
    det[i][2] = string of expansion, using standard 2,a,b,0 notation
    '''

    det = []
    coeffs = []
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
            # print('csfvec =', csfvec)
            # print('coeff = ', coeff)
            # print('a[i]  = ', a[i])
            # print('b[i]  = ', b[i])
            # print('i_a   = ', i_a)
            # print('i_b   = ', i_b)
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
            det.append(''.join(slat_det))
            coeffs.append(np.sqrt(np.abs(coeff)) * np.sign(coeff))

    return det, coeffs


def reading_molcas(file_to_read, root):
    CI = False
    EndCI = False
    closed=False
    CASSCF=True
    n_closed=1
    coeffs = []
    number = 0
    CIvec = []
    fin_CIvec = []
    fin_coeffs = []

    reading = False
    if CASSCF:
        with open(file_to_read) as f:
            for lines in f:

                if "printout of CI-coefficients larger than -1.00 for root  " + str(root) in lines:

                    reading = True
                elif "conf/sym" in lines and reading:
                    CI = True
                    EndCI = False

                elif "Natural orbitals and occupation numbers" in lines or " printout of CI-coefficients" in lines and CI:
                    EndCI = True
                    CI = False
                    break
                elif CI and not EndCI and lines != '\n' and lines:
                    s1 = lines.strip().split(" ")
                    s1= list(filter(None, s1))
                    print(s1)
                    CIvec.append(s1[1])
                    coeffs.append(float(s1[2]))
                    print(number)
                    number += 1


        coeffs = np.asarray(coeffs)
        print(sum(coeffs**2))
        print(coeffs)
        print(number)

        for i, j in zip(CIvec, coeffs):
            civec2, coeff2 = csf2sd(i, 0)
            coeff2 = np.asarray(coeff2)

            final_coeffs = j * coeff2
            for j, k in zip(civec2, final_coeffs):
                fin_CIvec.append('ab'*n_closed+j.replace('0','00').replace('a','a0').replace('b','0b').replace('2','ab'))

                fin_coeffs.append(k)

    else:

        with open(file_to_read) as f:
            for lines in f:

                if "The CI coefficients for the MIXED state nr.   " + str(root) in lines:
                    print(lines)
                    reading = True
                elif " Conf  SGUGA" in lines and reading:
                    CI = True
                    EndCI = False

                elif "The CI coefficients" in lines or "--" in lines and CI and reading:
                    EndCI = True
                    CI = False
                    break
                elif CI and not EndCI and lines != '\n' and lines:
                    s1 = lines.strip().split(" ")
                    s1 = list(filter(None, s1))
                    print(s1)
                    if len(s1)==8:
                        CIvec.append(s1[5])
                        coeffs.append(float(s1[6]))
                        print(number)
                        number += 1
                    else:
                        CIvec.append(s1[4])
                        coeffs.append(float(s1[5]))
                        print(number)
                        number += 1

        coeffs = np.asarray(coeffs)
        print(sum(coeffs ** 2))
        print(coeffs)
        print(number)
        number=0
        for i, j in zip(CIvec, coeffs):
            #print(number)
            #print(i)
            civec2, coeff2 = csf2sd(i, 0)
            coeff2 = np.asarray(coeff2)
            number+=1
            final_coeffs = j * coeff2
            for j, k in zip(civec2, final_coeffs):
                fin_CIvec.append(
                    'ab' * n_closed + j.replace('0', '00').replace('a', 'a0').replace('b', '0b').replace('2', 'ab'))

                fin_coeffs.append(k)
        print(fin_CIvec)
        print(fin_coeffs)

    return fin_CIvec, fin_coeffs

