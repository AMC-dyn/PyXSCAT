import numpy as np


def uunique(a):
    b, ipos, irep = np.unique(a, return_index=True, return_inverse=True, axis=1)
    print(np.size(ipos))
    sipos = np.sort(ipos)

    bsort = a[:, sipos]

    sirep = np.zeros(np.size(a[0, :]), dtype=np.int32)

    for i in range(np.size(b[0, :])):
        x = []
        # x = np.where((a == bsort[:, i]).all(axis=1))
        for j in range(np.size(a[0, :])):
            if np.all(a[:, j] == bsort[:, i]):
                x.append(j)
        x = np.asarray(x)
        sirep[x] = i


    sirep = np.asarray(sirep)
    sipos = np.asarray(sipos)
    return bsort, sipos, sirep
