import numpy as np

def uunique(a):
    b, ipos, irep = np.unique(a, return_index=True, return_inverse=True, axis=0)
    sipos = np.sort(ipos)
    bsort = a[sipos,:] 
    sirep = np.zeros(np.size(a[:,0]), dtype=np.int32)
    
    for i in range(np.size(b[:,0])):
        x = np.where((a == bsort[i,:]).all(axis=1))
        sirep[x] = i
    
    return bsort, sipos, sirep

