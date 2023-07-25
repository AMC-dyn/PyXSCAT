# -*- coding: utf-8 -*-
"""IAM generalized function

Automatically generated by Colaboratory.

Original file is located at
    https://colab.research.google.com/drive/1F4O6amCiNKKN6-y2R8hr95LYNfqOIYp0
"""

import numpy as np
import math
import matplotlib.pyplot as plt

# creating a dictionary to store data values in key:value pairs
elements = {
    "H": {"a": [0.493002, 0.322912, 0.140191, 0.040810],
          "b": [10.5109, 26.1257, 3.14236, 57.7997],
          "c": 0.003038
          },

    "Li": {"a": [1.12820, 0.750800, 0.617500, 0.465300],
           "b": [3.95460, 1.05240, 85.3905, 168.261],
           "c": 0.037700
           },

    "Li_plus": {"a": [0.696800, 0.788800, 0.341400, 0.156300],
                "b": [4.62370, 1.95570, 0.631600, 10.0953],
                "c": 0.016700
                },

    "C": {"a": [2.26069, 1.56165, 1.05075, 0.839259],
          "b": [22.6907, 0.656665, 9.75618, 55.5949],
          "c": 0.286977
          },

    "C_plus": {"a": [2.26069, 1.56165, 1.05075, 0.839259],
               "b": [22.6907, 0.656665, 9.75618, 55.5949],
               "c": 0.286977
               },

    "Ni": {"a": [12.2126, 3.13220, 2.01250, 1.16630],
           "b": [.005700, 9.89330, 28.9975, .582600],
           "c": -11.529
           },

    "F": {"a": [3.53920, 2.64120, 1.51700, 1.02430],
          "b": [10.2825, 4.29440, 0.261500, 26.1476],
          "c": 0.277600
          },

    "F_minus": {"a": [3.63220, 3.51057, 1.26064, 0.940706],
                "b": [5.27756, 14.7353, 0.442258, 47.3437],
                "c": 0.653396
                },

    "S": {
        "a": [6.90530, 5.20340, 1.43790, 1.58630],
        "b": [1.46790, 22.2151, 0.253600, 56.1720],
        "c": 0.866900},

    "O": {"a": [3.0485, 2.2868, 1.5463, 0.867],
          "b": [13.2771, 5.7011, 0.3239, 32.9089],
          "c": 0.2508
          },

    "I": {
        "a": [20.1472, 18.9949, 7.51380, 2.27350],
        "b": [4.34700, .381400, 27.7660, 66.8776],
        "c": 4.07120
    },

    "Xe": {
        "a": [20.2933, 19.0298, 8.97670, 1.99000],
        "b": [3.92820, 0.344000, 26.4659, 64.2658],
        "c": 3.71180}
}


def FF(a, b, c, q):
    # sum = return to 1x3 matrix by adding up the columns of 3x3 matrix
    expp = np.exp(-b * ((q / (4 * math.pi)) ** 2))
    array1 = a * expp
    return c + np.sum(array1, axis=0)


# geometry=molecular geometry,nAtom=number of atoms in the fragment,*args=name of the atom and number of that atom in the fragment such as ("H",2),("O",1),("F",4)
def plot_fnc(geometry, nAtoms, *args):
    counter = 0
    Nq = 100
    q = np.linspace(0.00000001, 7, num=100).reshape(1, 100)
    form_factor = np.zeros((nAtoms, Nq))

    # allocate intensity
    I = np.zeros((1, Nq))
    I_not = np.zeros((1, Nq))
    I_atom = np.zeros((1, Nq))

    for i in args:
        if i[0] in elements.keys():
            elm = i[0]
            FF_all = FF(np.asarray(elements[elm]["a"]).reshape(4, 1), np.asarray(elements[elm]["b"]).reshape(4, 1),
                        elements[elm]["c"], q)

            for j in range(i[1]):
                form_factor[counter, :] = np.asarray(FF_all)
                counter = +1

    # molecular scattering
    for i in range(1, nAtoms + 1):
        for j in range(1, nAtoms + 1):
            g = np.subtract(geometry[i - 1, :], geometry[j - 1, :])
            R_ij = np.sqrt(np.sum(np.power(g, 2), axis=0))

            s = np.transpose(np.squeeze(form_factor[i - 1, :]) * np.squeeze(form_factor[j - 1, :]))
            snc = np.sinc(q * R_ij / math.pi)
            I = I + np.multiply(s, snc)

    # create a figure object
    fig = plt.figure()

    y = np.transpose(I)
    x = np.transpose(q)
    print(y[0])
    # plot the functions
    plt.plot(x, y, 'c', label='IAM')
    plt.xlabel("q (1/Å)")
    plt.ylabel("I(q)")
    plt.legend(loc='upper right')

    # show the plot
    plt.show()


# geometry of the molecule
geometryH2 = np.array([[0.0, 0.0, 1.7],
                       [0.0, 0.0, 0.0]])

geometryHF = np.array([[0.0, 0.0, 0.917],
                       [0.0, 0.0, 0.0]])

plot_fnc(geometryH2, 2, ("H", 2))

plot_fnc(geometryHF, 2, ("O", 1), ("O", 1))
