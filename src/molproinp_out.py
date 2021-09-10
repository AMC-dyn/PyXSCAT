import os
import os.path
import re
import time

import numpy as np


def create_input():
    if os.path.isfile('../inputs/abinitio.dat'):
        print('reading abinitio input file')
    else:
        print('try to create your own input or change the name of the bloody file to the correct one, '
              'remember ----> abinitio.dat')
    RUN = False
    GEOM = False
    METHOD = False
    OCCUPATION = False
    CLOSED = False
    BASISlog = False
    NEL = False
    MLD = False
    PUN = False
    NUMSTATES = False
    MULTIPLICITY = False
    varname1 = 'n'
    varname2 = 'n'
    geomind = 0
    j = 0
    calctype = 'none'
    runningmolprocalc = 'def'
    with open('../inputs/abinitio.dat', 'r') as f:
        for lines in f:

            if lines.startswith('GEOMETRY'):
                GEOM = True
            if GEOM and not lines.startswith('GEOMETRY'):
                if geomind == 0:
                    units = lines.strip().split()
                    geomind = 1
                elif geomind == 1:
                    natoms = lines.strip().split()[0]
                    natoms = int(str(natoms))
                    atomname = []
                    x = np.zeros(natoms)
                    y = np.zeros(natoms)
                    z = np.zeros(natoms)
                    geomind = 2

                elif geomind == 2 and j <= natoms - 1:
                    atomname.append(lines.strip().split()[0])
                    x[j] = lines.strip().split()[1]
                    y[j] = lines.strip().split()[2]
                    z[j] = lines.strip().split()[3]
                    j = j + 1
                elif j > natoms - 1:
                    GEOM = False
            if lines.startswith('RUN'):
                RUN = True
            if not lines.startswith('RUN') and RUN:
                RUN = False
                runningmolprocalc = str(lines.strip().split()[0])

            if 'n' in runningmolprocalc or 'N' in runningmolprocalc:



                if lines.startswith('MLD'):
                    MLD = True

                if not lines.startswith('MLD') and MLD:
                    print('here')
                    MLD = False
                    varname1 = str(lines.strip().split()[0])
                    print(varname1)
                if lines.startswith('PUNCH'):
                    PUN = True
                if not lines.startswith('PUNCH') and PUN:
                    PUN = False
                    varname2 = str(lines.strip().split()[0])
                    print(varname2)

                if varname1 != 'n' and varname2 != 'n':
                    print('Copying the files')
                    # varname1=input('Specify the name of the molden file')
                    os.system('cp ' + varname1 + ' molpro.mld')
                    # varname2 = input('Specify the name of the punch')
                    os.system('cp ' + varname2 + ' molpro.pun')
                    return x, y, z
            elif 'y' in runningmolprocalc or 'Y' in runningmolprocalc:
                os.system('rm -f molpro.mld molpro.pun molpro*.out molpro*.xml ')
                if lines.startswith('METHOD'):
                    METHOD = True
                if not lines.startswith('METHOD') and METHOD:
                    calctype = str(lines.strip().split()[0])
                    METHOD = False
                if lines.startswith('OCC') and calctype == 'CAS' or calctype == 'cas':
                    OCCUPATION = True
                if lines.startswith('MULT'):
                    MULTIPLICITY = True
                if OCCUPATION and not lines.startswith('OCC'):
                    occ = str(lines.strip().split())[2]

                    occ = int(occ)
                    OCCUPATION = False
                if lines.startswith('CLOSED') and calctype == 'CAS' or calctype == 'cas':
                    CLOSED = True
                if CLOSED and not lines.startswith('CLOSED'):
                    closed = str(lines.strip().split())

                    closed = closed[2]
                    CLOSED = False
                if lines.startswith('BASIS'):
                    BASISlog = True
                if BASISlog and not lines.startswith('BASIS'):
                    basis = str(lines.strip().split(" "))
                    BASISlog = False
                if lines.startswith("NEL"):
                    NEL = True
                if NEL and not lines.startswith('NEL'):
                    NEL = False
                    numel = str(lines.strip())
                if lines.startswith('NSTATES'):
                    NUMSTATES = True
                if NUMSTATES and not lines.startswith('NSTATES'):
                    NUMSTATES = False
                    nstates = str(lines.strip())
                if MULTIPLICITY and not lines.startswith('MULT'):
                    MULTIPLICITY = False
                    if lines.strip() == 'singlet' or lines.strip() == 'Singlet' or lines.strip() == 'SINGLET':
                        mult = 0
                    elif lines.strip() == 'doublet' or lines.strip() == 'Doublet' or lines.strip() == 'DOUBLET':
                        mult = 1
                    else:
                        mult = 2

    # Molpro input, it must be changed for other codes
    '''routine to create molpro input, valid for molpro2012'''

    file = 'molpro_inp_scat_.inp'

    # if first:
    # if istep==0 and trajN==0:
    #     r=np.transpose(np.asarray([[7.0544201503E-01,-8.8768894715E-03,-6.6808940143E-03],[-6.8165975936E-01,1.7948934206E-02,4.0230273972E-03],[ 1.2640943417E+00,9.0767471618E-01,-3.0960211126E-02],[-1.4483835697E+00 ,8.7539956319E-01 ,4.6789959947E-02],[1.0430033100E+00,-9.0677165411E-01 ,8.1418967247E-02],[-1.1419988770E+00,-9.8436525752E-01,-6.5589218426E-02]]))
    with open(file, 'w') as f:
        f.write('***\n')
        f.write('memory,100,m\n')
        f.write('''gprint,civector,orbitals,angles=-1,distance=-1
 gthresh,twoint=1.0d-13
 gthresh,energy=1.0d-7,gradient=1.0d-2\n
 gthresh,thrpun=0.000000001
 ''')

        f.write('punch,molpro.pun,new\n')
        f.write('basis= ' + basis[2:-2] + '\n')

        f.write('''symmetry,nosym;
orient,noorient;
angstrom;
geomtype=xyz;
geom={
        ''')
        f.write(str(natoms) + '\n')
        f.write('\n')

        for i in range(natoms):
            line_wr = atomname[i] + ' ' + str(x[i]) + ' ' + str(y[i]) + ' ' + str(z[i]) + '\n'
            # print(file)
            # print(line_wr)
            f.write(line_wr)
        f.write('}\n')
        f.write('''{multi,failsafe;
maxiter,40;
config,det
''')
        line_wr = 'occ,' + str(occ) + ';\n'
        line_wr2 = 'closed,' + str(closed) + ';\n'
        line_wr3 = 'wf,' + str(numel) + ',1,' + str(mult) + ';\n'
        line_wr4 = 'state,' + str(nstates) + ';\n'

        f.write(line_wr)
        f.write(line_wr2)
        f.write(line_wr3)
        f.write(line_wr4)

        f.write('''pspace,10.0        
orbital,2140.3;
ORBITAL,IGNORE_ERROR;
ciguess,2501.2 
save,ci=2501.2}

put,molden,molpro.mld
''')

        f.write('''---''')

    os.system(
        'E:/Molpro/bin/molpro.exe -d . -s molpro_inp_scat_.inp')  # running molpro, change it for any run in a different computer
    time_counter = 0  #
    time_to_wait = 100
    while not os.path.exists(
            'molpro.pun'):  # Dodgy way to wait for the punch file to be created, must be other way more elegant
        time.sleep(1)
        time_counter += 1
        if time_counter > time_to_wait:
            print('more than 1000 secs')
            break
    print('Molpro input created, proceeding to read the variables from the molden and output files')
    return


# Starting the reading part, boring and annoying to program but very necessary....or not...;p
def outputreading(x, y, z):
    GTO = False
    MO = False
    SORB = False
    PORB = False
    DORB = False
    FORB = False
    GORB = False
    orbs = False
    realnum = []
    typ = []
    atoms = []
    rr = 1
    ga = []
    ci = []
    lang = []
    mang = []
    nang = []
    syms = []
    occs = []
    angpart = []
    MOS = []
    MOenergies = []
    MOnums = []
    xx = []
    yy = []
    zz = []
    cc = 0
    with open('molpro.mld') as f:
        for lines in f:

            if '[GTO]' in lines:
                GTO = True
            if '[MO]' in lines:
                MO = True
                GTO = False
            if GTO and not MO and not '[GTO]' in lines and len(lines) != 0 and len(lines.strip().split()) > 0:
                s1 = lines.strip().split(" ")

                s1 = [s for s in s1 if s]

                if '.' not in s1[0] and not re.fullmatch(r'^[a-z]+', s1[0]):
                    at = s1[0]
                elif re.fullmatch(r'^[a-z]+', s1[0]) is not None:

                    type = s1[0]
                    ncont = int(s1[1])

                    if type == 's':

                        orbs = True
                        SORB = True
                        ccc = 0
                    elif type == 'p':
                        orbs = True
                        PORB = True
                        ccc = 0

                    elif type == 'd':
                        orbs = True
                        DORB = True
                        ccc = 0
                    elif type == 'f':
                        orbs = True
                        FORB = True
                        ccc = 0
                    elif type == 'g':
                        orbs = True
                        GORB = True
                        ccc = 0
                    else:
                        print('This orbital type is not programmed yet')

                else:

                    if SORB and orbs:
                        xx.append(x[int(at) - 1])
                        yy.append(y[int(at) - 1])
                        zz.append(z[int(at) - 1])

                        realnum.append(cc + 1)
                        typ.append(type)
                        ga.append(float(s1[0].replace('D', 'E')))
                        ci.append(float(s1[1].replace('D', 'E')))
                        atoms.append(at)
                        lang.append(0)
                        mang.append(0)
                        nang.append(0)
                        angpart.append(rr)
                        ccc += 1
                        rr += 1
                        if ccc == ncont:
                            SORB = False
                            orbs = False
                            ccc = 0
                            cc += 1
                    if PORB and orbs:
                        xx.append(x[int(at) - 1])
                        yy.append(y[int(at) - 1])
                        zz.append(z[int(at) - 1])
                        realnum.append(cc + 1)
                        typ.append(type)
                        angpart.append(rr)
                        ga.append(float(s1[0].replace('D', 'E')))
                        ci.append(float(s1[1].replace('D', 'E')))
                        atoms.append(at)
                        lang.append(1)
                        mang.append(0)
                        nang.append(0)
                        xx.append(x[int(at) - 1])
                        yy.append(y[int(at) - 1])
                        zz.append(z[int(at) - 1])
                        realnum.append(cc + 2)
                        typ.append(type)
                        ga.append(float(s1[0].replace('D', 'E')))
                        ci.append(float(s1[1].replace('D', 'E')))
                        atoms.append(at)
                        angpart.append(rr)
                        lang.append(0)
                        mang.append(1)
                        nang.append(0)
                        xx.append(x[int(at) - 1])
                        yy.append(y[int(at) - 1])
                        zz.append(z[int(at) - 1])
                        realnum.append(cc + 3)
                        typ.append(type)
                        ga.append(float(s1[0].replace('D', 'E')))
                        ci.append(float(s1[1].replace('D', 'E')))
                        atoms.append(at)
                        lang.append(0)
                        mang.append(0)
                        nang.append(1)
                        angpart.append(rr)
                        ccc += 1
                        rr += 1
                        if ccc == ncont:
                            PORB = False
                            orbs = False
                            ccc = 0
                            cc += 3
                    elif DORB and orbs:
                        xx.append(x[int(at) - 1])
                        yy.append(y[int(at) - 1])
                        zz.append(z[int(at) - 1])
                        realnum.append(cc + 1)
                        typ.append(type)
                        angpart.append(rr)
                        ga.append(float(s1[0].replace('D', 'E')))
                        ci.append(float(s1[1].replace('D', 'E')))
                        atoms.append(at)
                        lang.append(2)
                        mang.append(0)
                        nang.append(0)

                        xx.append(x[int(at) - 1])
                        yy.append(y[int(at) - 1])
                        zz.append(z[int(at) - 1])
                        realnum.append(cc + 2)
                        typ.append(type)
                        ga.append(float(s1[0].replace('D', 'E')))
                        ci.append(float(s1[1].replace('D', 'E')))
                        atoms.append(at)
                        angpart.append(rr)
                        lang.append(0)
                        mang.append(2)
                        nang.append(0)

                        xx.append(x[int(at) - 1])
                        yy.append(y[int(at) - 1])
                        zz.append(z[int(at) - 1])
                        realnum.append(cc + 3)
                        typ.append(type)
                        ga.append(float(s1[0].replace('D', 'E')))
                        ci.append(float(s1[1].replace('D', 'E')))
                        atoms.append(at)
                        lang.append(0)
                        mang.append(0)
                        nang.append(2)
                        angpart.append(rr)

                        xx.append(x[int(at) - 1])
                        yy.append(y[int(at) - 1])
                        zz.append(z[int(at) - 1])
                        realnum.append(cc + 4)
                        typ.append(type)
                        ga.append(float(s1[0].replace('D', 'E')))
                        ci.append(float(s1[1].replace('D', 'E')))
                        atoms.append(at)
                        lang.append(1)
                        mang.append(1)
                        nang.append(0)
                        angpart.append(rr)

                        xx.append(x[int(at) - 1])
                        yy.append(y[int(at) - 1])
                        zz.append(z[int(at) - 1])
                        realnum.append(cc + 5)
                        typ.append(type)
                        ga.append(float(s1[0].replace('D', 'E')))
                        ci.append(float(s1[1].replace('D', 'E')))
                        atoms.append(at)
                        lang.append(1)
                        mang.append(0)
                        nang.append(1)
                        angpart.append(rr)

                        xx.append(x[int(at) - 1])
                        yy.append(y[int(at) - 1])
                        zz.append(z[int(at) - 1])
                        realnum.append(cc + 6)
                        typ.append(type)
                        ga.append(float(s1[0].replace('D', 'E')))
                        ci.append(float(s1[1].replace('D', 'E')))
                        atoms.append(at)
                        lang.append(0)
                        mang.append(1)
                        nang.append(1)
                        angpart.append(rr)

                        ccc += 1
                        rr += 1
                        if ccc == ncont:
                            DORB = False
                            orbs = False
                            ccc = 0
                            cc += 6

                    elif FORB and orbs:
                        xx.append(x[int(at) - 1])
                        yy.append(y[int(at) - 1])
                        zz.append(z[int(at) - 1])
                        realnum.append(cc + 1)
                        typ.append(type)
                        angpart.append(rr)
                        ga.append(float(s1[0].replace('D', 'E')))
                        ci.append(float(s1[1].replace('D', 'E')))
                        atoms.append(at)
                        lang.append(3)
                        mang.append(0)
                        nang.append(0)

                        xx.append(x[int(at) - 1])
                        yy.append(y[int(at) - 1])
                        zz.append(z[int(at) - 1])
                        realnum.append(cc + 2)
                        typ.append(type)
                        ga.append(float(s1[0].replace('D', 'E')))
                        ci.append(float(s1[1].replace('D', 'E')))
                        atoms.append(at)
                        angpart.append(rr)
                        lang.append(0)
                        mang.append(3)
                        nang.append(0)

                        xx.append(x[int(at) - 1])
                        yy.append(y[int(at) - 1])
                        zz.append(z[int(at) - 1])
                        realnum.append(cc + 3)
                        typ.append(type)
                        ga.append(float(s1[0].replace('D', 'E')))
                        ci.append(float(s1[1].replace('D', 'E')))
                        atoms.append(at)
                        lang.append(0)
                        mang.append(0)
                        nang.append(3)
                        angpart.append(rr)

                        xx.append(x[int(at) - 1])
                        yy.append(y[int(at) - 1])
                        zz.append(z[int(at) - 1])
                        realnum.append(cc + 4)
                        typ.append(type)
                        ga.append(float(s1[0].replace('D', 'E')))
                        ci.append(float(s1[1].replace('D', 'E')))
                        atoms.append(at)
                        lang.append(1)
                        mang.append(2)
                        nang.append(0)
                        angpart.append(rr)

                        xx.append(x[int(at) - 1])
                        yy.append(y[int(at) - 1])
                        zz.append(z[int(at) - 1])
                        realnum.append(cc + 5)
                        typ.append(type)
                        ga.append(float(s1[0].replace('D', 'E')))
                        ci.append(float(s1[1].replace('D', 'E')))
                        atoms.append(at)
                        lang.append(2)
                        mang.append(1)
                        nang.append(0)
                        angpart.append(rr)

                        xx.append(x[int(at) - 1])
                        yy.append(y[int(at) - 1])
                        zz.append(z[int(at) - 1])
                        realnum.append(cc + 6)
                        typ.append(type)
                        ga.append(float(s1[0].replace('D', 'E')))
                        ci.append(float(s1[1].replace('D', 'E')))
                        atoms.append(at)
                        lang.append(2)
                        mang.append(0)
                        nang.append(1)
                        angpart.append(rr)

                        xx.append(x[int(at) - 1])
                        yy.append(y[int(at) - 1])
                        zz.append(z[int(at) - 1])
                        realnum.append(cc + 7)
                        typ.append(type)
                        ga.append(float(s1[0].replace('D', 'E')))
                        ci.append(float(s1[1].replace('D', 'E')))
                        atoms.append(at)
                        lang.append(1)
                        mang.append(0)
                        nang.append(2)
                        angpart.append(rr)

                        xx.append(x[int(at) - 1])
                        yy.append(y[int(at) - 1])
                        zz.append(z[int(at) - 1])
                        realnum.append(cc + 8)
                        typ.append(type)
                        ga.append(float(s1[0].replace('D', 'E')))
                        ci.append(float(s1[1].replace('D', 'E')))
                        atoms.append(at)
                        lang.append(0)
                        mang.append(1)
                        nang.append(2)
                        angpart.append(rr)

                        xx.append(x[int(at) - 1])
                        yy.append(y[int(at) - 1])
                        zz.append(z[int(at) - 1])
                        realnum.append(cc + 9)
                        typ.append(type)
                        ga.append(float(s1[0].replace('D', 'E')))
                        ci.append(float(s1[1].replace('D', 'E')))
                        atoms.append(at)
                        lang.append(0)
                        mang.append(2)
                        nang.append(1)
                        angpart.append(rr)

                        xx.append(x[int(at) - 1])
                        yy.append(y[int(at) - 1])
                        zz.append(z[int(at) - 1])
                        realnum.append(cc + 10)
                        typ.append(type)
                        ga.append(float(s1[0].replace('D', 'E')))
                        ci.append(float(s1[1].replace('D', 'E')))
                        atoms.append(at)
                        lang.append(1)
                        mang.append(1)
                        nang.append(1)
                        angpart.append(rr)

                        ccc += 1
                        rr += 1
                        if ccc == ncont:
                            FORB = False
                            orbs = False
                            ccc = 0
                            cc += 10

                    elif GORB and orbs:
                        xx.append(x[int(at) - 1])
                        yy.append(y[int(at) - 1])
                        zz.append(z[int(at) - 1])
                        realnum.append(cc + 1)
                        typ.append(type)
                        angpart.append(rr)
                        ga.append(float(s1[0].replace('D', 'E')))
                        ci.append(float(s1[1].replace('D', 'E')))
                        atoms.append(at)
                        lang.append(4)
                        mang.append(0)
                        nang.append(0)

                        xx.append(x[int(at) - 1])
                        yy.append(y[int(at) - 1])
                        zz.append(z[int(at) - 1])
                        realnum.append(cc + 2)
                        typ.append(type)
                        ga.append(float(s1[0].replace('D', 'E')))
                        ci.append(float(s1[1].replace('D', 'E')))
                        atoms.append(at)
                        angpart.append(rr)
                        lang.append(0)
                        mang.append(4)
                        nang.append(0)

                        xx.append(x[int(at) - 1])
                        yy.append(y[int(at) - 1])
                        zz.append(z[int(at) - 1])
                        realnum.append(cc + 3)
                        typ.append(type)
                        ga.append(float(s1[0].replace('D', 'E')))
                        ci.append(float(s1[1].replace('D', 'E')))
                        atoms.append(at)
                        lang.append(0)
                        mang.append(0)
                        nang.append(4)
                        angpart.append(rr)

                        xx.append(x[int(at) - 1])
                        yy.append(y[int(at) - 1])
                        zz.append(z[int(at) - 1])
                        realnum.append(cc + 4)
                        typ.append(type)
                        ga.append(float(s1[0].replace('D', 'E')))
                        ci.append(float(s1[1].replace('D', 'E')))
                        atoms.append(at)
                        lang.append(3)
                        mang.append(1)
                        nang.append(0)
                        angpart.append(rr)

                        xx.append(x[int(at) - 1])
                        yy.append(y[int(at) - 1])
                        zz.append(z[int(at) - 1])
                        realnum.append(cc + 5)
                        typ.append(type)
                        ga.append(float(s1[0].replace('D', 'E')))
                        ci.append(float(s1[1].replace('D', 'E')))
                        atoms.append(at)
                        lang.append(3)
                        mang.append(0)
                        nang.append(1)
                        angpart.append(rr)

                        xx.append(x[int(at) - 1])
                        yy.append(y[int(at) - 1])
                        zz.append(z[int(at) - 1])
                        realnum.append(cc + 6)
                        typ.append(type)
                        ga.append(float(s1[0].replace('D', 'E')))
                        ci.append(float(s1[1].replace('D', 'E')))
                        atoms.append(at)
                        lang.append(1)
                        mang.append(3)
                        nang.append(0)
                        angpart.append(rr)

                        xx.append(x[int(at) - 1])
                        yy.append(y[int(at) - 1])
                        zz.append(z[int(at) - 1])
                        realnum.append(cc + 7)
                        typ.append(type)
                        ga.append(float(s1[0].replace('D', 'E')))
                        ci.append(float(s1[1].replace('D', 'E')))
                        atoms.append(at)
                        lang.append(0)
                        mang.append(3)
                        nang.append(1)
                        angpart.append(rr)

                        xx.append(x[int(at) - 1])
                        yy.append(y[int(at) - 1])
                        zz.append(z[int(at) - 1])
                        realnum.append(cc + 8)
                        typ.append(type)
                        ga.append(float(s1[0].replace('D', 'E')))
                        ci.append(float(s1[1].replace('D', 'E')))
                        atoms.append(at)
                        lang.append(1)
                        mang.append(0)
                        nang.append(3)
                        angpart.append(rr)

                        xx.append(x[int(at) - 1])
                        yy.append(y[int(at) - 1])
                        zz.append(z[int(at) - 1])
                        realnum.append(cc + 9)
                        typ.append(type)
                        ga.append(float(s1[0].replace('D', 'E')))
                        ci.append(float(s1[1].replace('D', 'E')))
                        atoms.append(at)
                        lang.append(0)
                        mang.append(1)
                        nang.append(3)
                        angpart.append(rr)

                        xx.append(x[int(at) - 1])
                        yy.append(y[int(at) - 1])
                        zz.append(z[int(at) - 1])
                        realnum.append(cc + 10)
                        typ.append(type)
                        ga.append(float(s1[0].replace('D', 'E')))
                        ci.append(float(s1[1].replace('D', 'E')))
                        atoms.append(at)
                        lang.append(2)
                        mang.append(2)
                        nang.append(0)
                        angpart.append(rr)

                        xx.append(x[int(at) - 1])
                        yy.append(y[int(at) - 1])
                        zz.append(z[int(at) - 1])
                        realnum.append(cc + 11)
                        typ.append(type)
                        ga.append(float(s1[0].replace('D', 'E')))
                        ci.append(float(s1[1].replace('D', 'E')))
                        atoms.append(at)
                        lang.append(2)
                        mang.append(0)
                        nang.append(2)
                        angpart.append(rr)

                        xx.append(x[int(at) - 1])
                        yy.append(y[int(at) - 1])
                        zz.append(z[int(at) - 1])
                        realnum.append(cc + 12)
                        typ.append(type)
                        ga.append(float(s1[0].replace('D', 'E')))
                        ci.append(float(s1[1].replace('D', 'E')))
                        atoms.append(at)
                        lang.append(0)
                        mang.append(2)
                        nang.append(2)
                        angpart.append(rr)

                        xx.append(x[int(at) - 1])
                        yy.append(y[int(at) - 1])
                        zz.append(z[int(at) - 1])
                        realnum.append(cc + 13)
                        typ.append(type)
                        ga.append(float(s1[0].replace('D', 'E')))
                        ci.append(float(s1[1].replace('D', 'E')))
                        atoms.append(at)
                        lang.append(2)
                        mang.append(1)
                        nang.append(1)
                        angpart.append(rr)

                        xx.append(x[int(at) - 1])
                        yy.append(y[int(at) - 1])
                        zz.append(z[int(at) - 1])
                        realnum.append(cc + 14)
                        typ.append(type)
                        ga.append(float(s1[0].replace('D', 'E')))
                        ci.append(float(s1[1].replace('D', 'E')))
                        atoms.append(at)
                        lang.append(1)
                        mang.append(2)
                        nang.append(1)
                        angpart.append(rr)

                        xx.append(x[int(at) - 1])
                        yy.append(y[int(at) - 1])
                        zz.append(z[int(at) - 1])
                        realnum.append(cc + 15)
                        typ.append(type)
                        ga.append(float(s1[0].replace('D', 'E')))
                        ci.append(float(s1[1].replace('D', 'E')))
                        atoms.append(at)
                        lang.append(1)
                        mang.append(1)
                        nang.append(2)
                        angpart.append(rr)

                        ccc += 1
                        rr += 1
                        if ccc == ncont:
                            GORB = False
                            orbs = False
                            ccc = 0
                            cc += 15

            if MO and '[MO]' not in lines and len(lines) != 0 and len(lines.strip().split()) > 0:

                if 'Sym' in lines:
                    s1 = lines.strip().split(" ")
                    s1 = [s for s in s1 if s]
                    syms.append(s1[1])
                elif 'Occ' in lines:
                    s1 = lines.strip().split(" ")
                    s1 = [s for s in s1 if s]
                    occs.append(float(s1[1]))
                elif 'Ene' in lines:
                    s1 = lines.strip().split(" ")
                    s1 = [s for s in s1 if s]
                    MOenergies.append(s1[1])
                elif 'Spin' not in lines:
                    s1 = lines.strip().split(" ")
                    s1 = [s for s in s1 if s]
                    MOS.append(s1[1])
                    MOnums.append(s1[0])
        # Number of active orbitals, in a calculation with SA the occupation of the virtual orbitals can be 0 or
        # close to 0

        total = len(occs)
        xx = np.asarray(xx, dtype=np.float64)
        yy = np.asarray(yy, dtype=np.float64)
        zz = np.asarray(zz, dtype=np.float64)

        ga = np.asarray(ga, dtype=np.float64)
        ci = np.asarray(ci, dtype=np.float64)

        lang = np.asarray(lang, dtype=np.int64)
        mang = np.asarray(mang, dtype=np.int64)
        nang = np.asarray(nang, dtype=np.int64)

        # It is necessary to reorder the orbitals to appear in the same order as they do in the output
        # We can use the symmetries to do this

        idx = np.argsort(np.asarray(syms, dtype=np.float64))

        mos = np.reshape(MOS, (total, max(realnum)))

        mos = mos[idx, :]
        occs = np.asarray(occs)[idx]
        actives = sum(map(lambda x: x != 0.0000, occs))

        # We restrict the number of orbitals to those with an occupation bigger than 0
        mos = mos[:actives, :]

        # Realnum just duplicates the MOs for the different contractions and/or the different angular momenta in the
        # basis sets
        realnum = np.asarray(realnum)
        mos = mos[:, realnum - 1]

    return ga, ci, realnum, mang, lang, nang, mos, MOnums, actives, total, angpart, xx, yy, zz
