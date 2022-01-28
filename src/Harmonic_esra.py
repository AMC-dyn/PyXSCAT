import numpy as np
import math
def fnc(t,k,m,x,v):
    w=np.sqrt(k/m)
    return x*math.cos(w*t)+(v/w)*math.sin(w*t)

time_vector=np.linspace(0, 15, num=100000)

x_position=[]
k=3
for t in time_vector:
  x_position.append((fnc(t,k,0.7,0,2)))
pos=np.asarray(x_position)

import matplotlib.pyplot as plt
import numpy as np

x = pos
y = k*x**2/2

# setting the axes at the centre
fig = plt.figure()
ax = fig.add_subplot(1, 1, 1)
ax.spines['left'].set_position('center')
ax.spines['bottom'].set_position('center')
ax.spines['right'].set_color('none')
ax.spines['top'].set_color('none')
ax.xaxis.set_ticks_position('bottom')
ax.yaxis.set_ticks_position('left')
# plot the functions
plt.plot(x,time_vector, 'c', label='f')


plt.legend(loc='upper left')

# show the plot
plt.show()
