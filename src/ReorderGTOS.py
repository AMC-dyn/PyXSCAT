import numpy as np


def reorder(l, m, n, ga, ci, M, angpart):
    apos, apos2, arep = np.unique(angpart, return_index=True, return_inverse=True)


    N_new = len(apos)
    counter = 0

    ord_vec = np.zeros(len(l))


    for i in range(N_new):
        duplicates_ang = np.argwhere(arep == arep[apos2[i]])
        # print(duplicates_ang)
        num = len(duplicates_ang)


        if num == 1:
            ord_vec[counter] = duplicates_ang[0]
            counter += 1

        elif num == 3:
            ord_vec[counter] = duplicates_ang[2]
            ord_vec[counter + 1] = duplicates_ang[1]
            ord_vec[counter + 2] = duplicates_ang[0]

            counter += 3

        elif num == 6:
            ord_vec[counter] = duplicates_ang[2]
            ord_vec[counter + 1] = duplicates_ang[5]
            ord_vec[counter + 2] = duplicates_ang[1]
            ord_vec[counter + 3] = duplicates_ang[4]
            ord_vec[counter + 4] = duplicates_ang[3]
            ord_vec[counter + 5] = duplicates_ang[0]
            counter = counter + 6
        elif num == 10:
            ord_vec[counter] = duplicates_ang[2]
            ord_vec[counter + 1] = duplicates_ang[7]
            ord_vec[counter + 2] = duplicates_ang[8]
            ord_vec[counter + 3] = duplicates_ang[1]
            ord_vec[counter + 4] = duplicates_ang[6]
            ord_vec[counter + 5] = duplicates_ang[9]
            ord_vec[counter + 6] = duplicates_ang[3]
            ord_vec[counter + 7] = duplicates_ang[5]
            ord_vec[counter + 8] = duplicates_ang[4]
            ord_vec[counter + 9] = duplicates_ang[0]
            counter = counter + 10
        elif num == 15:
            ord_vec[counter] = duplicates_ang[2]
            ord_vec[counter + 1] = duplicates_ang[7]
            ord_vec[counter + 2] = duplicates_ang[8]
            ord_vec[counter + 3] = duplicates_ang[1]
            ord_vec[counter + 4] = duplicates_ang[6]
            ord_vec[counter + 5] = duplicates_ang[9]
            ord_vec[counter + 6] = duplicates_ang[3]
            ord_vec[counter + 7] = duplicates_ang[5]
            ord_vec[counter + 8] = duplicates_ang[4]
            ord_vec[counter + 9] = duplicates_ang[0]
            ord_vec[counter + 10] = duplicates_ang[9]
            ord_vec[counter + 11] = duplicates_ang[3]
            ord_vec[counter + 12] = duplicates_ang[5]
            ord_vec[counter + 13] = duplicates_ang[4]
            ord_vec[counter + 14] = duplicates_ang[0]

    ord_vec = np.array(ord_vec, dtype=np.int16)

    l = np.asarray(l)
    m = np.asarray(m)
    n = np.asarray(n)
    ga = np.asarray(ga)
    ci = np.asarray(ci)
    M = np.asarray(M, dtype=np.float64)
    angpart = np.asarray(angpart)



    return l[ord_vec], m[ord_vec], n[ord_vec], ga[ord_vec], ci[ord_vec], M[:, ord_vec], angpart[
        ord_vec]
