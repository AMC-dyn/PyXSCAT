import sys
from textwrap import wrap

import numpy as np

import molcas_ci_reader as mc
import molden_reader_nikola as mldreader
import molden_reader_nikola_bagel as bgmldreader
import molden_reader_nikola_molcas as mcmldreader
import twordm_red as td

# If there is an external file with CIvecs or 2rdms JeremyR==True

input_dict = {}
input_file = sys.argv[1]

try:
    with open(input_file, 'r') as f:
        for line in f:
            if line[0] != '#':
                if '#' in line:
                    line = line.split('#')[0]
                key, value = line.strip('\n').strip(' ').split('=', 1)
                input_dict[key.upper()] = value
except FileNotFoundError:
    print('No input file, reverting to defaults!')

try:
    state1 = int(input_dict['STATE1'])
except KeyError:
    state1 = 1
    print(f'No state specified, performing scattering on state {state1}')

try:
    state2 = int(input_dict['STATE2'])
except KeyError:
    state2 = state1

try:
    closed = int(input_dict['CLOSED'])
    print(f'Number of closed orbitals = {closed}')
except ValueError:
    closed = [int(i) for i in input_dict['CLOSED'].split()]
except KeyError:
    print('Not given "closed", assuming no closed orbitals')
    closed = 0

try:
    qmin, qmax, npoints = input_dict['QVEC'].split()
    qmin = float(qmin)
    if qmin < 1E-10:
        qmin = 1E-10
        print(f"rounding up qmin to {qmin} for numerical reasons")
    qmax = float(qmax)
    npoints = int(npoints)
except KeyError:
    print("No q vector specified, defaulting to")
    qmin = 1E-10
    qmax = 10
    npoints = 100

print(f"{npoints} points between {qmin} AA and {qmax} AA")

try:
    input_dict['TYPE']
except KeyError:
    print('No "TYPE" provided, defaulting to total scattering')
    input_dict['TYPE'] = 'TOTAL'

try:
    input_dict['PARTICLE']
except KeyError:
    print('No "PARTICLE" provided, defaulting to x-ray scattering')
    input_dict['PARTICLE'] = 'XRAY'

try:
    input_dict['BESSELORDER']
except KeyError:
    print('No "BESSELORDER" provided, defaulting to J0 scattering')
    input_dict['BESSELORDER'] = '0'

input_dict['TYPE'] = input_dict['TYPE'].upper()
input_dict['PARTICLE'] = input_dict['PARTICLE'].upper()

if input_dict['TYPE'] == 'TOTAL':
    if input_dict['PARTICLE'] == 'XRAY':
        if input_dict['BESSELORDER'] == '0':
            Type = 1
        elif input_dict['BESSELORDER'] == '2':
            Type = 7
    elif input_dict['PARTICLE'] == 'ELECTRON':
        if input_dict['BESSELORDER'] == '0':
            Type = 5
        else:
            print(
                'Electron scattering is currently not available except with J0 scattering'
            )
            sys.exit()
elif input_dict['TYPE'] == 'ELASTIC':
    if input_dict['PARTICLE'] == 'XRAY':
        if input_dict['BESSELORDER'] == '0':
            Type = 2
        elif input_dict['BESSELORDER'] == '2':
            Type = 8
    elif input_dict['PARTICLE'] == 'ELECTRON':
        if input_dict['BESSELORDER'] == '0':
            Type = 6
        else:
            print(
                'Electron scattering is currently not available except with J0 scattering'
            )
            sys.exit()
elif input_dict['TYPE'] == 'IAM':
    print('Will soon be implemented!!!')
    sys.exit()

if state2 == state1:
    title = f"J{input_dict['BESSELORDER']} {input_dict['TYPE']} {input_dict['PARTICLE']} scattering calculation for state {state1}"
elif input_dict['TYPE'] == 'TOTAL':
    title = f"J{input_dict['BESSELORDER']} {input_dict['PARTICLE']} coherent mixed scattering calculation between states {state1} and {state2}"
elif input_dict['TYPE'] == 'ELASTIC':
    title = f"J{input_dict['BESSELORDER']} {input_dict['PARTICLE']} inelastic scattering calculation between states {state1} and {state2}"
else:
    print("Not sure here... exiting...")
    sys.exit()

print(title.capitalize())

try:
    mldfile = input_dict['MOFILE']
except KeyError:
    print('You must provide a file with molecular orbitals, exiting...')
    sys.exit()

try:
    method = input_dict['METHOD'].upper()
except KeyError:
    print('No method provided, assuming Hartree-Fock')

if method == 'CI':
    ci = True
    hf = False
    try:
        punfile = input_dict['CIFILE']
    except KeyError:
        print('You must provide a file with a ci vector, exiting...')
        sys.exit()
    print(
        f"Using CI style expansion, reading MOs from {mldfile} and CI vector from {punfile}"
    )
else:
    ci = False
    hf = True
    print(f"Using HF style wavefunction, reading MOs from {mldfile}")

try:
    outfile = input_dict['OUTFILE']
except KeyError:
    outfile = 'scattering.out'

print(f"Printing scattering data to {outfile}")

# TODO
# Currently just using defaults, will make input for these soon!!

cutoffcentre = 0.1  # suggested: cutoffcentre = 0.01;
# the cutoff for the Z integral
cutoffz = 1e-20  # suggested: cutoffz = 1E-9;
# the cutoff for the product of the MD coefficients
cutoffmd = 1e-20  # suggested: cutoffmd = 1E-20;
# Largeconf is used when the det space is very large but we have still no external file so we need to create one

largeconf = False
# Type of calculation
# TOTAL--> 1
# ELASTIC --> 2
# TOTAL ALIGNED --> 3
# ELASTIC ALIGNED --> 4
# TOTAL ELECTRON--> 5
# ELASTIC ELECTRON --> 6
# TOTAL J2 --> 7
# ELASTIC J2 --> 8
# Ouput name

if ci:
    # read ci vectors
    print(f"Reading CI vectors from {punfile}")
    civs, confs = td.twordmconst(closed, punfile)

    Nmo_max = len(confs[0]) / 2
    print(f"{Nmo_max} occupied orbitals")

elif hf:
    # create single slater determinant form
    print(f"This needs to be fixed, currently only 5 orbitals")
    civs = [1.000]
    # The number of occupied orbs or he configuration used must be specified by the user in hf
    norbs = 5
    confs = ['ab' * norbs]
    print(confs)
    Nmo_max = 100
else:
    print(f"Not sure what this does...")
    confs = 0
    civs = 0
    Nmo_max = 100
civs = np.array(civs)

#read molden file

print(f"Reading {mldfile}")
Nmo_max = 100
gtos, atoms = mldreader.read_orbitals(mldfile, N=Nmo_max, decontract=True)
geom = atoms.geometry()


#Write output

with open('options.dat', 'w') as f:
    f.write(str(np.size(atoms.atomic_numbers())) + '\n')
    for i in atoms.atomic_numbers():
        f.write(str(i) + ' ')
    f.write('\n')
    for i in range(np.size(atoms.atomic_numbers())):
        f.write(str(geom[i, :])[1:-1] + '\n')
    f.write(str(cutoffcentre) + '\n')
    f.write(str(cutoffz) + '\n')
    f.write(str(cutoffmd) + '\n')
    f.write(str(False) + '\n')
    f.write(str(False) + '\n')
    f.write(str(hf) + '\n')
    f.write(str(qmin) + ' ' + str(qmax) + ' ' + str(npoints) + ' \n')
    f.write(str(Type) + '\n')
    f.write(str(state1) + ' ' + str(state2) + '\n')
    f.write(outfile + '\n')
    f.write(str(True) + '\n')
    f.write(str(False) + '\n')
    f.write(str(False) + '\n')
    if ci:
        f.write(str(np.size(confs)) + '\n')

        f.write(str(len(confs[0])) + '\n')
        f.write(str(np.size(civs[0, :])) + '\n')

        if np.size(confs) < 1E4:
            bitwise = False
            print('Normal integration')
            for i in range(np.size(confs)):

                # [*confs[i].replace('a', '1').replace('b', '2')]
                lst = wrap(confs[i].replace('a', '1').replace('b', '2'), 1)
                confs2 = np.asarray(lst, dtype=np.int64)

                g1 = [int(i) for i in lst]
                f.write(str(g1)[1:-1].replace(',', ' '))
                for j in range(np.size(civs[0, :])):
                    f.write(' ' + str(civs[i, j]) + ' ')
                f.write(str('\n'))
        else:
            print('bitwise integration')
            f.write('bitwise.dat')
            bitwise = True
            alpha = 0
            beta = 0
            with open('bitwise.dat', 'w') as m:
                for i in range(np.size(confs)):

                    lst = [*confs[i].replace('a', '1').replace('b', '2')]
                    confs2 = np.asarray(lst, dtype=np.int64)
                    alpha = 0
                    beta = 0
                    for j in range(np.size(confs2)):
                        if confs2[j] != 0:
                            if np.mod(j, 2) == 0:
                                alpha = alpha + 2**(j / 2)
                            else:
                                beta = beta + 2**((j - 1) / 2)
                    m.write(
                        str(i) + ' ' + str(civs[i, state1 - 1]) + ' ' +
                        str(civs[i, state2 - 1]) + ' ' + str(int(alpha)) +
                        ' ' + str(int(beta)) + '\n')
        f.write(str(bitwise) + '\n')
    elif hf:
        f.write(str(np.size(confs)) + '\n')

        f.write(str(len(confs[0])) + '\n')

        f.write(str(1) + '\n')

        lst = [*confs[0].replace('a', '1').replace('b', '2')]
        confs2 = np.asarray(lst, dtype=np.int64)
        f.write(str(confs2)[1:-1])
        f.write(' ' + str(1.00000) + ' ')
        f.write(str('\n'))

Nmo_max = 100
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
mmod = np.asarray(mmod, dtype=np.float64)

with open('basis.dat', 'w') as f:
    f.write(str(np.size(l)) + '\n')
    for i in range(np.size(l)):
        f.write(
            str(xx[i]) + ' ' + str(yy[i]) + ' ' + str(zz[i]) + ' ' +
            str(ga[i]) + ' ' + str(l[i]) + ' ' + str(m[i]) + ' ' + str(n[i]) +
            ' ' + str(group[i]) + '\n')
with open('MOs.dat', 'w') as f:
    f.write(str(np.size(mmod[:, 0])) + ' ' + str(np.size(mmod[0, :])) + '\n')
    for i in range(np.size(mmod[:, 0])):
        for j in range(np.size(mmod[0, :])):
            f.write(str(mmod[i, j]) + ' ')
        f.write('\n')
