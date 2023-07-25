"""Define molecules and their atoms.

Author: Nikola Zotev
Last modified: 30/03/2020
"""

import numpy as np
import pdb

AU2ANG = 0.52917721092


class Atom:
    """Symbol, number, coordinates and mass of atoms."""

    def __init__(self, atmnum, x, y, z, xyz_units='au', charge=0):
        """Construct an atom.

        symbol -- str -- atomic symbol
        atmnum -- int -- atomic number
        x, y, z -- float --- coordinates
        charge  - int   --- the charge, used to set ions
        """
        self.atmnum = atmnum
        self.x = x
        self.y = y
        self.z = z
        # by default it is a neutral atom not ion
        self.charge = charge
        # connected to atom type
        self.symbol = self.atom_symbol(atmnum)
        # self.mass = selatmnumf.atomic_mass(atmnum)
        self.xyz_units = xyz_units

    def set_default_mass(self):
        """Use to set the mass of the atom the the default value in amu."""
        self.mass = self.atomic_mass(self.atmnum)
        self.mass_units = 'amu'

    def set_mass(self, mass, mass_units='amu'):
        """Use to set the mass if different from the default."""
        self.mass = mass
        self.mass_units = mass_units

    @staticmethod
    def atom_symbol(atmnum):
        """Set symbol of atom."""
        return symbol_dict[atmnum]

    @staticmethod
    def atomic_mass(atmnum):
        """Set atomic mass in amu.

        Which isotopes to use as default?
        """
        return mass_dict[atmnum]

    def mass_in_kg(self):
        """Mass in kilograms."""
        return self.mass * 1.66053904E-27

    def mass_in_au(self):
        """Mass in au of mass."""
        return self.mass * 1822.888486209


class Molecule:
    """Define molecule as a collection of atoms."""

    def __init__(self):
        """Make an empty molecule."""
        self.atoms = []
        self.n = 0

    def add_atom(self, new_atom):
        """Add atom to the molecule."""
        # modify to make sure that all atoms use the same units of distance
        self.atoms.append(new_atom)
        self.n += 1

    def center_of_mass(self):
        """Find the cener of mass."""
        pass

    def move_to_com(self):
        """Move to the the center of mass."""
        pass

    def atomic_numbers(self):
        """Return atomic numbers as a numpy array."""
        atmnums = []
        for i in self.atoms:
            atmnums.append(i.atmnum)

        return np.array(atmnums)

    def geometry(self):
        """Return the geometry as a 2D np array."""
        geom = np.zeros((self.n, 3))
        for i, at in enumerate(self.atoms):
            geom[i, :] = [at.x, at.y, at.z]

        return geom


def xyz_geometry_read(xyzfile):
    """Read molecule from *.xyz file."""
    with open(xyzfile) as file:
        # the first line contain the number of atoms
        line = file.readline()
        Nat = int(line)
        # skip the next line - it text message
        line = file.readline()
        # or molecule
        mol = Molecule()
        for i in range(0, Nat):
            line = file.readline()
            lsplit = line.split()
            symbol = lsplit[0]
            # as tupid way to use a dictionarry backwards
            for a, sym in symbol_dict.items():
                if sym == symbol:
                    break
            x = float(lsplit[1])
            y = float(lsplit[2])
            z = float(lsplit[3])
            atom = Atom(a, x, y, z, xyz_units='ang')
            mol.add_atom(atom)

    return mol

mass_dict = {1  : 1.008 ,
           2  : 4.003 ,
           3  : 6.941 ,
           4  : 9.012 ,
           5  : 10.811,
           6  : 12.011,
           7  : 14.007,
           8  : 15.999,
           9  : 18.998,
           10 : 20.180,
           11 : 22.990,
           12 : 24.305,
           13 : 26.982,
           14 : 28.086,
           15 : 30.974,
           16 : 32.065,
           17 : 35.453,
           18 : 39.948,
           19 : 39.098,
           20 : 40.078,
           21 : 44.956,
           22 : 47.867,
           23 : 50.942,
           24 : 51.996,
           25 : 54.938,
           26 : 55.845,
           27 : 58.933,
           28 : 58.693,
           29 : 63.546,
           30 : 65.390,
           31 : 69.723,
           32 : 72.640,
           33 : 74.922,
           34 : 78.960,
           35 : 79.904,
           36 : 83.800,
           37 : 85.468,
           38 : 87.620,
           39 : 88.906,
           40 : 91.224,
           41 : 92.906,
           42 : 95.940,
           43 : 98.000,
           44 : 101.070,
           45 : 102.906,
           46 : 106.420,
           47 : 107.868,
           48 : 112.411,
           49 : 114.818,
           50 : 118.710,
           51 : 121.760,
           52 : 127.600,
           53 : 126.905,
           54 : 131.293,
           55 : 132.906,
           56 : 137.327,
           57 : 138.906,
           58 : 140.116,
           59 : 140.908,
           60 : 144.240,
           61 : 145.000,
           62 : 150.360,
           63 : 151.964,
           64 : 157.250,
           65 : 158.925,
           66 : 162.500,
           67 : 164.930,
           68 : 167.259,
           69 : 168.934,
           70 : 173.040,
           71 : 174.967,
           72 : 178.490,
           73 : 180.948,
           74 : 183.840,
           75 : 186.207,
           76 : 190.230,
           77 : 192.217,
           78 : 195.078,
           79 : 196.967,
           80 : 200.590,
           81 : 204.383,
           82 : 207.200,
           83 : 208.980,
           84 : 209.000,
           85 : 210.000,
           86 : 222.000,
           87 : 223.000,
           88 : 226.000,
           89 : 227.000,
           90 : 232.038,
           91 : 231.036,
           92 : 238.029,
           93 : 237.000,
           94 : 244.000,
           95 : 243.000,
           96 : 247.000,
           97 : 247.000,
           98 : 251.000,
           99 : 252.000,
           100: 257.000,
           101: 258.000,
           102: 259.000,
           103: 262.000,
           104: 261.000,
           105: 262.000,
           106: 266.000,
           107: 264.000,
           108: 277.000,
           109: 268.000
           } # end mass_dict

symbol_dict = {
                1  : "H",
                2  : "He",
                3  : "Li",
                4  : "Be",
                5  : "B",
                6  : "C",
                7  : "N",
                8  : "O",
                9  : "F",
                10 : "Ne",
                11 : "Na",
                12 : "Mg",
                13 : "Al",
                14 : "Si",
                15 : "P",
                16 : "S",
                17 : "Cl",
                18 : "Ar",
                19 : "K",
                20 : "Ca",
                21 : "Sc",
                22 : "Ti",
                23 : "V",
                24 : "Cr",
                25 : "Mn",
                26 : "Fe",
                27 : "Co",
                28 : "Ni",
                29 : "Cu",
                30 : "Zn",
                31 : "Ga",
                32 : "Ge",
                33 : "As",
                34 : "Se",
                35 : "Br",
                36 : "Kr",
                37 : "Rb",
                38 : "Sr",
                39 : "Y",
                40 : "Zr",
                41 : "Nb",
                42 : "Mo",
                43 : "Tc",
                44 : "Ru",
                45 : "Rh",
                46 : "Pd",
                47 : "Ag",
                48 : "Cd",
                49 : "In",
                50 : "Sn",
                51 : "Sb",
                52 : "Te",
                53 : "I",
                54 : "Xe",
                55 : "Cs",
                56 : "Ba",
                57 : "La",
                58 : "Ce",
                59 : "Pr",
                60 : "Nd",
                61 : "Pm",
                62 : "Sm",
                63 : "Eu",
                64 : "Gd",
                65 : "Tb",
                66 : "Dy",
                67 : "Ho",
                68 : "Er",
                69 : "Tm",
                70 : "Yb",
                71 : "Lu",
                72 : "Hf",
                73 : "Ta",
                74 : "W",
                75 : "Re",
                76 : "Os",
                77 : "Ir",
                78 : "Pt",
                79 : "Au",
                80 : "Hg",
                81 : "Tl",
                82 : "Pb",
                83 : "Bi",
                84 : "Po",
                85 : "At",
                86 : "Rn",
                87 : "Fr",
                88 : "Ra",
                89 : "Ac",
                90 : "Th",
                91 : "Pa",
                92 : "U",
                93 : "Np",
                94 : "Pu",
                95 : "Am",
                96 : "Cm",
                97 : "Bk",
                98 : "Cf",
                99 : "Es",
                100: "Fm",
                101: "Md",
                102: "No",
                103: "Lr",
                104: "Rf",
                105: "0",
                106: "Sg",
                107: "0",
                108: "0",
                109: "Mt"
} # end dict_sybmol
