"""Reading components of molden files.

This is more of a objected oriented design but one can probably do this with
structured numpy arrays - more similar to to MATLAB/Fortran style. The syntax
of structured arrays is not intuitive to me so prefer to encapsulate the same
functionality in objects.
"""

import numpy as np
import molecule as mol
import re
import math
import copy
import pdb

AU2ANG = 0.52917721092


class PrimitiveGTO:
    """Describe primitive GTO."""

    def __init__(self, atom_idx, contraction, group, ga, c, l, m, n):
        """Define Primitive.

        atom_idx - int - the atom on which the GTO is centered
        ga - float - the exponent
        c - float - the contraction coefficient
        contraction - int - which contraction does it belong to
        l, m, n - int - the angular momentum q. numbers
        group - int - e.g. px, py and pz component of a p GTO form a group
                    (useful for MD coefficients etc)
        """
        self.atom_idx = atom_idx
        self.ga = ga
        self.c = c
        self.contraction = contraction
        self.l = l
        self.m = m
        self.n = n
        self.group = group
        # MO coefficients
        self.mo = np.empty(1)
        self.x = None
        self.y = None
        self.z = None

    def print_gto(self):
        """Print GTO."""
        temp = (self.atom_idx, self.contraction, self.group, self.ga, self.c,
                self.l, self.m, self.n)


    def norm(self):
        """Normalise GTO."""
        normA = (2 / math.pi) ** 0.75
        normB = 2 ** (self.l + self.m + self.n)
        normC = self.ga ** ((2 * self.l + 2 * self.m + 2 * self.n + 3) / 4.0)
        normD = factd(2 * self.l - 1) * factd(2 * self.m - 1) * factd(2 * self.n - 1)

        norm = normA * normB * normC / normD ** 0.5
        return norm


class ArrayGTO:
    """Describe orbitals and GTOs."""

    def __init__(self, mo, x, y, z, ga, l, m, n, group):
        """Read all arguments as numpy arrays."""
        self.mo = mo
        self.x = x
        self.y = y
        self.z = z
        self.ga = ga
        self.l = l
        self.m = m
        self.n = n
        self.group = group


def _read_atoms_file(file):
    """Read the atomic coordinates and numbers from mld file.

    input:
    file -- open file to read from
    output:
    atoms -- list of type Atom
    """

    # start from the beginning of the file
    file.seek(0)
    # find the posisiton of the tag [Atoms]
    for line in file:
        if "[Atoms]" in line:
            break
    else:
        print("Something went wrong! Can't find [Atoms].")

    # reading the atomic positions
    atoms = mol.Molecule()

    for line in file:
        if "[GTO]" in line:
            break
        else:
            columns = line.split()
            # symbol = columns[0]
            atmnum = int(columns[2])
            # coordinates in au
            x = float(columns[3])/AU2ANG
            y = float(columns[4])/AU2ANG
            z = float(columns[5])/AU2ANG
            atoms.add_atom(mol.Atom(atmnum, x, y, z))

    else:
        print("Something went wrong! Can't find [GTO]")

    return atoms


def read_geometry(file):
    """Read molecule from molden file.

    This is wrapper. Can take both file object or filename.
    Returns an instance of class Molecule.
    """
    if isinstance(file, str):
        with open(file, 'r') as f:
            return _read_atoms_file(f)
    else:
        return _read_atoms_file(file)


def _read_contractions(file):
    """Read the contracted GTOs.

    Reads the GTOs, and expands the angular momentum part
    input:
        file -- open file
    output:
        gtos -- list of instances of PrimitiveGTO
    """
    # start from the beginning of the file
    file.seek(0)
    # find the the GTO tag [GTO]
    for line in file:
        if "[GTO]" in line:
            break
    else:
        print("Something went wrong! Can't find [GTO].")
        return

    # starting writing the GTOs into lists
    GTOs = []

    # read the first atom
    line = file.readline()
    atms = [int(s) for s in line.split() if s.isdigit()]

    atm = atms[0]



    contraction_counter = 1
    group_counter = 0
    # loop until you find the [MO] tag
    for line in file:
        if "[MO]" in line:
            break
        else:
            # skip if line is blank or whitspace
            if not line or line.isspace():

                line = file.readline()
                # if the second line is blank or whitespace just go the the MO
                # condition
                if not line or line.isspace():
                    continue
                # if it was a single blank line read the next atom
                else:
                    atms = [int(s) for s in line.split() if s.isdigit()]
                    atm = atms[0]

                    line = file.readline()
            # read the contraction
            contraction_spec = line.split()

            type_of_GTO = contraction_spec[0]
            num_of_primitives = int(contraction_spec[1])

            for i in range(0, num_of_primitives):
                line = file.readline()
                line = line.replace("D", "E")
                g = float(line.split()[0])
                c = float(line.split()[1])
                group_counter += 1

                # dealing with the angular momentum
                ang = []
                if type_of_GTO == 's':
                    ang.append([0, 0, 0])
                elif type_of_GTO == 'p':
                    ang.append([1, 0, 0])
                    ang.append([0, 1, 0])
                    ang.append([0, 0, 1])
                elif type_of_GTO == 'd':
                    ang.append([2, 0, 0])
                    ang.append([0, 2, 0])
                    ang.append([0, 0, 2])
                    ang.append([1, 1, 0])
                    ang.append([1, 0, 1])
                    ang.append([0, 1, 1])
                elif type_of_GTO == 'f':
                    ang.append([3, 0, 0])
                    ang.append([0, 3, 0])
                    ang.append([0, 0, 3])
                    ang.append([1, 2, 0])
                    ang.append([2, 1, 0])
                    ang.append([2, 0, 1])
                    ang.append([1, 0, 2])
                    ang.append([0, 1, 2])
                    ang.append([0, 2, 1])
                    ang.append([1, 1, 1])
                else:
                    print("Error! I can only read s, p, d and f GTOs.")
                # all ang momentum for the given contraction
                for i, a in enumerate(ang):
                    l = a[0]
                    m = a[1]
                    n = a[2]
                    temp_gto = PrimitiveGTO(atm, contraction_counter + i,
                                            group_counter, g, c, l, m, n)
                    GTOs.append(temp_gto)
                    temp_gto.print_gto()

            contraction_counter += len(ang)
    else:
        print("Something went wrong! Can't find [MO]")
    print('size GTOS', np.size(GTOs))
    return GTOs


def _read_MO(file, mo_cutoff):
    """Read the MOs from molden file.

    input:
        file -- open file
        mo_cutoff -- how many orbitals to _read_atoms_file
    output:
        mo - numpy matrix that contains the MO coeffs
        occ_array - numpy array that contains the occupancy of the orbitals
        energy_array - numpy array - the energy of the orbitals
        syms_array - numpy array - the symmetry labels of the orbitals
    """
    # start from the beginning of the file
    file.seek(0)
    # find the the GTO tag [GTO]
    for line in file:
        if "[MO]" in line:
            break
    else:
        print("Something went wrong! Can't find [MO].")
        return
    # symmetry labels
    syms = []
    # energy
    energy = []
    # spin
    spin = []
    # occupancies = useful for single state natural orbitals
    occ = []
    # the MO coefficients
    mo_table = []
    mo_counter = 0
    mo_this = []

    # reading the MO, first energies and occupancies than the coeffs
    line = file.readline()
    while line and mo_counter <= mo_cutoff and "[FREQ]" not in line:
        if "Sym" in line:
            mo_counter += 1
            if mo_counter <= mo_cutoff:
                s = re.findall("\d+\.\d+", line)
                if s:
                    syms.append(s[0])
                else:
                    syms.append(1.1)
            if mo_this:
                mo_table.append(mo_this)
            # the coefficient for this orbitals
            mo_this = []
        elif "Ene" in line:
            s = re.findall("\d+\.\d+", line)
            energy.append(float(s[0]))
        elif "Spin" in line:
            s = line.split()
            spin.append(s[1].strip())
        elif "Occup" in line:
            oc = re.findall("\d+\.\d+", line)
            occ.append(float(oc[0]))
        else:
            temp = line.split()
            mo_this.append(float(temp[1]))
        # move to a new line after you record the current one
        line = file.readline()

    # record the last orbital before hitting the end of the file
    if mo_this:
        mo_table.append(mo_this)

    # redored the MO according to the Molpro labels
    nMO = len(syms)
    syms_array = np.array([float(i) for i in syms])
    idx1 = np.argsort(syms_array)
    print('previous order', syms)
    syms_array = (syms_array - syms_array.astype(int)) * 1000 + syms_array.astype(int)

    idx = np.argsort(syms_array)


    mo = np.array(mo_table)
    syms_array = np.array([float(i) for i in syms])
    print('after order', syms_array[idx])

    mo = mo[idx, :]

    energy_array = np.array(energy)
    energy_array = energy_array[idx]

    occ_array = np.array(occ)
    occ_array = occ_array[idx]
    print('final occupation is', sum(occ_array), 'with ', np.size(idx))

    return (np.transpose(mo), occ_array, energy_array, syms_array)


def factd(n):
    """Do double factorial."""
    if n <= 0:
        return 1
    return n * factd(n - 2)


def _mo_fill_gto(GTOs, mo_table):
    """Add the MO for each GTOs.

    input:
        GTOs - list of PrimitiveGTO - table of GTOs
        mo_table - numpy matrix - mo coefficients
    """
    for gto in GTOs:
        gto.mo = mo_table[gto.contraction - 1, :]



def _reoder_gto():
    """Reorder GTOs angular momentum."""
    # do I need that? This is needed inside the integral in a very specific
    # order. Maybe there is a better solution though
    pass


def _reverse_contraction_scheme(GTOs):
    """Reduce basis set by removing contractions."""
    # The list of contracted GTOs
    contracted = []
    # add empty GTO to start
    temp_gto = copy.deepcopy(GTOs[0])
    temp_gto.mo = 0
    contracted.append(temp_gto)
    for i in GTOs:
        for j in contracted:
            # if the primitive is already included as a part of another
            # contraction
            if (i.atom_idx == j.atom_idx and i.ga == j.ga
                    and i.l == j.l and i.m == j.m and i.n == j.n):
                # add the MO to the total (assumed that the gto has been
                # normalised)
                j.mo += i.mo
                break
        # if it does not find it anywhere
        else:
            this_gto = copy.deepcopy(i)
            # Remove the obsolted attributes
            this_gto.contraction = None
            this_gto.c = None
            # add the current one
            this_gto.mo = i.mo
            contracted.append(this_gto)

    # redo the numbering of the groups so that it goes 1,2,3... without breaks
    counter = 1
    prev_group = contracted[0].group

    for i in range(1, len(contracted)):
        if (contracted[i].group != prev_group):
            counter += 1
            prev_group = contracted[i].group
            contracted[i].group = counter
        else:
            contracted[i].group = counter

    return contracted


def _normalise_gto(GTOs):
    """Normalise the GTO so that the MO coeffs are become weigths."""
    for i in GTOs:
        i.mo = i.c * i.norm() * i.mo
        i.c = None


def _xyz_fill_gto(GTOs, molecule):
    """Fill the position of the atoms."""
    for i in GTOs:
        i.x = molecule.atoms[i.atom_idx - 1].x
        i.y = molecule.atoms[i.atom_idx - 1].y
        i.z = molecule.atoms[i.atom_idx - 1].z


def _gtos2numpy(GTOs):
    """Return numpy tables and arrays for the GTOs."""
    nGTOs = len(GTOs)
    nMO = GTOs[0].mo.size

    l = np.zeros(nGTOs, dtype=int, order='F')
    m = np.zeros(nGTOs, dtype=int, order='F')
    n = np.zeros(nGTOs, dtype=int, order='F')
    ga = np.zeros(nGTOs, dtype=float, order='F')
    x = np.zeros(nGTOs, dtype=float, order='F')
    y = np.zeros(nGTOs, dtype=float, order='F')
    z = np.zeros(nGTOs, dtype=float, order='F')
    group = np.zeros(nGTOs, dtype=int, order='F')

    mo = np.zeros((nGTOs, nMO), dtype=float, order='F')

    for i in range(0, nGTOs):
        l[i] = GTOs[i].l
        m[i] = GTOs[i].m
        n[i] = GTOs[i].n
        ga[i] = GTOs[i].ga
        x[i] = GTOs[i].x
        y[i] = GTOs[i].y
        z[i] = GTOs[i].z
        group[i] = GTOs[i].group
        mo[i, :] = GTOs[i].mo

    return ArrayGTO(mo, x, y, z, ga, l, m, n, group)


def _read_MO_GTOs(file, N, decontract=False):
    """Read the MO and GTOs from an open file."""
    # read the geometry
    molecule = _read_atoms_file(file)
    # read the GTOs
    GTOs = _read_contractions(file)
    # read the MO and assign to each primitive
    (mo_table, occ, energy, syms) = _read_MO(file, N)

    _mo_fill_gto(GTOs, mo_table)
    _normalise_gto(GTOs)

    # reverse the contraction scheme to leave only unrepeated primitives
    if decontract:
        GTOs = _reverse_contraction_scheme(GTOs)
    # assign the xyz of the atoms
    _xyz_fill_gto(GTOs, molecule)
    # return numpy arrays
    gto_array = _gtos2numpy(GTOs)
    gto_array.occ = occ
    gto_array.energy = energy
    gto_array.syms = syms

    return gto_array


def read_orbitals(file, N=10000, decontract=False):
    """Read the GTOs and MOs.

    This is wrapper. Can take both file object or filename.
    Returns an tuple of numpy arrays.
    """
    if isinstance(file, str):
        with open(file, 'r') as f:
            gtos = _read_MO_GTOs(f, N, decontract)
            atoms = _read_atoms_file(f)
    else:
        return _read_MO_GTOs(file, N, decontract)
    return gtos, atoms


def read_occupied_orbitals(file, decontract=False):
    """Read only the occupied orbitals."""
    N = 10000
    if isinstance(file, str):
        with open(file, 'r') as f:
            gto_array = _read_MO_GTOs(f, N, decontract)
    else:
        gto_array = _read_MO_GTOs(file, N, decontract)

    slicer = np.where(gto_array.occ == 0)[0]
    gto_array.occ = np.delete(gto_array.occ, slicer)
    gto_array.energy = np.delete(gto_array.energy, slicer)
    gto_array.mo = np.delete(gto_array.mo, slicer, 1)
    gto_array.syms = np.delete(gto_array.syms, slicer)

    return gto_array

# if __name__ == "__main__":
#     with open('./INPUT/NH3/nh3_phf_avdz.mld', 'r') as f:
#         atoms = _read_atoms_file(f)
#         GTOs = _read_contractions(f)
#         (mo_table, occ, energy) = _read_MO(f, 5)
#         _mo_fill_gto(GTOs, mo_table)
#         gtos_pure = _reverse_contraction_scheme(GTOs)
#         _xyz_fill_gto(gtos_pure, atoms)
#         _gtos2numpy(gtos_pure)
#
