import numpy as np


def reorder(l, m, n, ga, ci, M, angpart):
    apos, apos2, arep = np.unique(angpart, return_index=True, return_inverse=True)
    print(arep)
    print(apos - 1)
    N_new = len(apos)
    counter = 0
    print(N_new)
    ord_vec = np.zeros(len(l))
    print(len(l))
    print(angpart)
    for i in range(N_new):
        duplicates_ang = np.argwhere(arep == arep[apos2[i]])
        # print(duplicates_ang)
        num = len(duplicates_ang)
        print(arep[apos2[i]], num)

        if num == 1:
            ord_vec[counter] = duplicates_ang[0]
            counter += 1

        elif num == 3:
            ord_vec[counter] = duplicates_ang[2]
            ord_vec[counter + 1] = duplicates_ang[1]
            ord_vec[counter + 2] = duplicates_ang[0]

            counter += 3

        elif num == 6:
            ord_vec[counter] = duplicates_ang[5]
            ord_vec[counter + 1] = duplicates_ang[4]
            ord_vec[counter + 2] = duplicates_ang[3]
            ord_vec[counter + 3] = duplicates_ang[2]
            ord_vec[counter + 4] = duplicates_ang[1]
            ord_vec[counter + 5] = duplicates_ang[0]
            counter = counter + 6
    print('thats the ord_vec', ord_vec)
    ord_vec = np.array(ord_vec, dtype=np.int16)

    l = np.asarray(l)
    m = np.asarray(m)
    n = np.asarray(n)
    ga = np.asarray(ga)
    ci = np.asarray(ci)
    M = np.asarray(M, dtype=np.float32)
    angpart = np.asarray(angpart)

    print(l[ord_vec])

    return l[ord_vec], m[ord_vec], n[ord_vec], ga[ord_vec], ci[ord_vec], M[:, ord_vec], angpart[
        ord_vec]
