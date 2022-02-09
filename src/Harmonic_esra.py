
import numpy as np
import matplotlib.pyplot as plt

# general function for form factor (f):
# f = @(a,b,c,q) c + sum(a .* exp(-b *(q/(4*pi))^2))
# a: vector, b: vector, c: scalar, q: scattering vector
a_H_t = np.array([0.493002, 0.322912, 0.140191,  0.040810]).reshape(4,1)
b_H_t = np.array([10.5109,  26.1257,   3.14236,  57.7997]).reshape(4,1)
c_H = 0.003038
q=np.linspace(0.00000001,7,num=100)




import math
def FF(a,b,c,q):
  #sum = Conversion to 1x3 matrix by adding up the columns of 3x3 matrix.
  expp=np.exp(-b*((q/(4*math.pi))**2))
  array1=a*expp
  return c+np.sum(array1,axis=0)

q_new=np.array([[0.000000000001, 0.07070707, 0.14141414, 0.21212121, 0.28282828,
       0.35353535, 0.42424242, 0.49494949, 0.56565657, 0.63636364,
       0.70707071, 0.77777778, 0.84848485, 0.91919192, 0.98989899,
       1.06060606, 1.13131313, 1.2020202 , 1.27272727, 1.34343434,
       1.41414141, 1.48484848, 1.55555556, 1.62626263, 1.6969697 ,
       1.76767677, 1.83838384, 1.90909091, 1.97979798, 2.05050505,
       2.12121212, 2.19191919, 2.26262626, 2.33333333, 2.4040404 ,
       2.47474747, 2.54545455, 2.61616162, 2.68686869, 2.75757576,
       2.82828283, 2.8989899 , 2.96969697, 3.04040404, 3.11111111,
       3.18181818, 3.25252525, 3.32323232, 3.39393939, 3.46464646,
       3.53535354, 3.60606061, 3.67676768, 3.74747475, 3.81818182,
       3.88888889, 3.95959596, 4.03030303, 4.1010101 , 4.17171717,
       4.24242424, 4.31313131, 4.38383838, 4.45454545, 4.52525253,
       4.5959596 , 4.66666667, 4.73737374, 4.80808081, 4.87878788,
       4.94949495, 5.02020202, 5.09090909, 5.16161616, 5.23232323,
       5.3030303 , 5.37373737, 5.44444444, 5.51515152, 5.58585859,
       5.65656566, 5.72727273, 5.7979798 , 5.86868687, 5.93939394,
       6.01010101, 6.08080808, 6.15151515, 6.22222222, 6.29292929,
       6.36363636, 6.43434343, 6.50505051, 6.57575758, 6.64646465,
       6.71717172, 6.78787879, 6.85858586, 6.92929293, 7.        ]]) #(1,100)


#You need to define the form-factors as a matrix of dimension (natoms x qpoints). In this case you have 2 atoms only-> form_factor needs to have dimension 2 x 100
#Try to avoid using the same name for functions and variables, it can create a mess in the code

form_factor=np.zeros((2,100))
for i in range(0,2):
    form_factor[i,:]=np.asarray(FF(a_H_t,b_H_t,c_H,q_new))

print(form_factor[0])
#geometry of H2
geometry=np.array([[0.0 , 0.0 , 1.7],
                   [0.0 , 0.0 , 0.0]])

Nq=100 #length of q space
nAtoms=2 #for H2 (number of atoms in the fragment)

I=np.zeros((1,Nq))
I_not=np.zeros((1,Nq))
I_atom =np.zeros((1,Nq))

for i in range(1,nAtoms+1):
  for j in range(1,nAtoms+1):
    g=np.subtract(geometry[i-1,:],geometry[j-1,:])
    R_ij=np.sqrt(np.sum(np.power(g,2),axis=0))

    s=np.transpose(np.squeeze(form_factor[i-1,:])*np.squeeze(form_factor[j-1,:]))
    snc=np.sinc(q_new*R_ij/math.pi)
    I= I + np.multiply(s,snc)




# setting the axes at the centre
fig = plt.figure()

y=np.transpose(I)
x=np.transpose(q_new)

# plot the functions
plt.plot(x,y, 'c', label='f')
plt.xlabel("q (1/Å)")
plt.ylabel("I")
plt.legend(loc='upper right')

# show the plot
plt.show()