import numpy as np

atnum = 2
geom = np.zeros((atnum, 3))

geom[0, 2] = 1.7

q = np.linspace(0, 10, 100)

x = 0
for i in range(10):
    for j in range(2):
        x = x + 1
        if i==8:
            print(i)



