import numpy as np


def md_table_gen(l, x, ga):
    N = len(l)
    maxl = np.max(l)
    D = np.zeros((N, N, 2 * maxl + 1, maxl + 1, maxl + 1))
    print(2*maxl+1)

    for i in range(0, N):
        for j in range(i + 1, N):
            gaP = ga[i] + ga[j]
            Px = (ga[i] * x[i] + ga[j] * x[j]) / gaP

            if l[i] < l[j]:
                ii = j
                jj = i
            else:
                ii = i
                jj = j

            PA = Px - x[ii]
            PB = Px - x[jj]

            l1 = l[ii]
            l2 = l[jj]
            if l1 == 0 and l2 == 0:
                D[ii, jj, 0, 0, 0] = 1
                D[jj, ii, 0, 0, 0] = 1
            elif l1 == 1 and l2 == 0:
                D[ii, jj, 0, 1, 0] = PA
                D[ii, jj, 1, 1, 0] = 0.5 / gaP
                D[jj, ii, 0, 0, 1] = D[ii, jj, 0, 1, 0]
                D[jj, ii, 1, 0, 1] = D[ii, jj, 1, 1, 0]

            elif l1 == 1 and l2 == 1:
                a = 0.5 / gaP
                D[ii, jj, 0, 0, 1] = PB
                D[ii, jj, 1, 0, 1] = a
                D[ii, jj, 0, 1, 0] = PA
                D[ii, jj, 1, 1, 0] = a

                D[ii, jj, 0, 1, 1] = PB * D[ii, jj, 0, 1, 0] + D[ii, jj, 1, 1, 0]
                D[ii, jj, 1, 1, 1] = a * D[ii, jj, 0, 1, 0] + PB * D[ii, jj, 1, 1, 0]
                D[ii, jj, 2, 1, 1] = a * D[ii, jj, 1, 1, 0]

                D[jj, ii, 0, 1, 1] = D[ii, jj, 0, 1, 1]
                D[jj, ii, 1, 1, 1] = D[ii, jj, 1, 1, 1]
                D[jj, ii, 2, 1, 1] = D[ii, jj, 2, 1, 1]
            elif l1 == 2 and l2 == 0:
                a = 0.5 / gaP

                D[ii, jj, 0, 1, 0] = PA
                D[ii, jj, 1, 1, 0] = a

                D[ii, jj, 0, 2, 0] = PA * D[ii, jj, 0, 1, 0] + D[ii, jj, 1, 1, 0]  # corrected
                D[ii, jj, 1, 2, 0] = a * D[ii, jj, 0, 1, 0] + PA * D[ii, jj, 1, 1, 0]  # corrected
                D[ii, jj, 2, 2, 0] = a * D[ii, jj, 1, 1, 0]  # corrected

                D[jj, ii, 0, 0, 2] = D[ii, jj, 0, 2, 0]
                D[jj, ii, 1, 0, 2] = D[ii, jj, 1, 2, 0]
                D[jj, ii, 2, 0, 2] = D[ii, jj, 2, 2, 0]
            elif l1 == 2 and l2 == 1:
                a = 0.5 / gaP
                D[ii, jj, 0, 0, 1] = PB
                D[ii, jj, 1, 0, 1] = a
                D[ii, jj, 0, 1, 0] = PA
                D[ii, jj, 1, 1, 0] = a
                D[ii, jj, 0, 1, 1] = PB * D[ii, jj, 0, 1, 0] + D[ii, jj, 1, 1, 0]  # corrected
                D[ii, jj, 1, 1, 1] = a * D[ii, jj, 0, 1, 0] + PB * D[ii, jj, 1, 1, 0]  # corrected
                D[ii, jj, 2, 1, 1] = a * D[ii, jj, 1, 1, 0]
                D[ii, jj, 0, 2, 0] = PA * D[ii, jj, 0, 1, 0] + D[ii, jj, 1, 1, 0]
                D[ii, jj, 1, 2, 0] = a * D[ii, jj, 0, 1, 0] + PA * D[ii, jj, 1, 1, 0]
                D[ii, jj, 2, 2, 0] = a * D[ii, jj, 1, 1, 0]

                D[ii, jj, 0, 2, 1] = PB * D[ii, jj, 0, 2, 0] + D[ii, jj, 1, 2, 0]
                D[ii, jj, 1, 2, 1] = a * D[ii, jj, 0, 2, 0] + PB * D[ii, jj, 1, 2, 0] + 2. * D[ii, jj, 2, 2, 0]
                D[ii, jj, 2, 2, 1] = a * D[ii, jj, 1, 2, 0] + PB * D[ii, jj, 2, 2, 0]
                D[ii, jj, 3, 2, 1] = a * D[ii, jj, 2, 2, 0]

                D[jj, ii, 0, 1, 2] = D[ii, jj, 0, 2, 1]
                D[jj, ii, 1, 1, 2] = D[ii, jj, 1, 2, 1]
                D[jj, ii, 2, 1, 2] = D[ii, jj, 2, 2, 1]
                D[jj, ii, 3, 1, 2] = D[ii, jj, 3, 2, 1]
            elif l1 == 2 and l2 == 2:
                a = 0.5 / gaP

                D[ii, jj, 0, 0, 1] = PB
                D[ii, jj, 1, 0, 1] = a
                D[ii, jj, 0, 0, 2] = PB * D[ii, jj, 0, 0, 1] + D[ii, jj, 1, 0, 1]
                D[ii, jj, 1, 0, 2] = a * D[ii, jj, 0, 0, 1] + PB * D[ii, jj, 1, 0, 1]
                D[ii, jj, 2, 0, 2] = a * D[ii, jj, 1, 0, 1]
                D[ii, jj, 0, 1, 0] = PA
                D[ii, jj, 1, 1, 0] = a  # corrected
                D[ii, jj, 0, 1, 1] = PB * D[ii, jj, 0, 1, 0] + D[ii, jj, 1, 1, 0]  # corrected
                D[ii, jj, 1, 1, 1] = a * D[ii, jj, 0, 1, 0] + PB * D[ii, jj, 1, 1, 0]  # corrected
                D[ii, jj, 2, 1, 1] = a * D[ii, jj, 1, 1, 0]  # corrected
                D[ii, jj, 0, 1, 2] = PB * D[ii, jj, 0, 1, 1] + D[ii, jj, 1, 1, 1]  # corrected
                D[ii, jj, 1, 1, 2] = a * D[ii, jj, 0, 1, 1] + PB * D[ii, jj, 1, 1, 1] + 2. * D[ii, jj, 2, 1, 1]  # corrected
                D[ii, jj, 2, 1, 2] = a * D[ii, jj, 1, 1, 1] + PB * D[ii, jj, 2, 1, 1]  # corrected
                D[ii, jj, 3, 1, 2] = a * D[ii, jj, 2, 1, 1]

                D[ii, jj, 0, 2, 0] = PA * D[ii, jj, 0, 1, 0] + D[ii, jj, 1, 1, 0]  # corrected
                D[ii, jj, 1, 2, 0] = a * D[ii, jj, 0, 1, 0] + PA * D[ii, jj, 1, 1, 0]  # corrected
                D[ii, jj, 2, 2, 0] = a * D[ii, jj, 1, 1, 0]  # corrected
                D[ii, jj, 0, 2, 1] = PB * D[ii, jj, 0, 2, 0] + D[ii, jj, 1, 2, 0]
                D[ii, jj, 1, 2, 1] = a * D[ii, jj, 0, 2, 0] + PB * D[ii, jj, 1, 2, 0] + 2. * D[ii, jj, 2, 2, 0]
                D[ii, jj, 2, 2, 1] = a * D[ii, jj, 1, 2, 0] + PB * D[ii, jj, 2, 2, 0]
                D[ii, jj, 3, 2, 1] = a * D[ii, jj, 2, 2, 0]
                D[ii, jj, 0, 2, 2] = PB * D[ii, jj, 0, 2, 1] + D[ii, jj, 1, 2, 1]
                D[ii, jj, 1, 2, 2] = a * D[ii, jj, 0, 2, 1] + PB * D[ii, jj, 1, 2, 1] + 2. * D[ii, jj, 2, 2, 1]
                D[ii, jj, 2, 2, 2] = a * D[ii, jj, 1, 2, 1] + PB * D[ii, jj, 2, 2, 1] + 3. * D[ii, jj, 3, 2, 1]
                D[ii, jj, 3, 2, 2] = a * D[ii, jj, 2, 2, 1] + PB * D[ii, jj, 3, 2, 1]
                D[ii, jj, 4, 2, 2] = a * D[ii, jj, 3, 2, 1]

                D[jj, ii, 0, 2, 2] = D[ii, jj, 0, 2, 2]
                D[jj, ii, 1, 2, 2] = D[ii, jj, 1, 2, 2]
                D[jj, ii, 2, 2, 2] = D[ii, jj, 2, 2, 2]
                D[jj, ii, 3, 2, 2] = D[ii, jj, 3, 2, 2]
                D[jj, ii, 4, 2, 2] = D[ii, jj, 4, 2, 2]

            elif l1 == 3 and l2 == 0:
                a = 0.5 / gaP
                D[ii, jj, 0, 1, 0] = PA
                D[ii, jj, 1, 2, 0] = a
                D[ii, jj, 0, 2, 0] = PA * D[ii, jj, 0, 1, 0] + D[ii, jj, 1, 2, 0]
                D[ii, jj, 1, 2, 0] = a * D[ii, jj, 0, 1, 0] + PA * D[ii, jj, 1, 2, 0]
                D[ii, jj, 2, 2, 0] = a * D[ii, jj, 1, 2, 0]
                D[ii, jj, 0, 3, 0] = PA * D[ii, jj, 0, 2, 0] + D[ii, jj, 1, 2, 0]
                D[ii, jj, 1, 3, 0] = a * D[ii, jj, 0, 2, 0] + PA * D[ii, jj, 1, 2, 0] + 2. * D[ii, jj, 2, 2, 0]
                D[ii, jj, 2, 3, 0] = a * D[ii, jj, 1, 2, 0] + PA * D[ii, jj, 2, 2, 0]
                D[ii, jj, 3, 3, 0] = a * D[ii, jj, 2, 2, 0]

                D[jj, ii, 0, 0, 3] = D[ii, jj, 0, 3, 0]
                D[jj, ii, 1, 0, 3] = D[ii, jj, 1, 3, 0]
                D[jj, ii, 2, 0, 3] = D[ii, jj, 2, 3, 0]
                D[jj, ii, 3, 0, 3] = D[ii, jj, 3, 3, 0]

            elif l1 == 3 and l2 == 1:
                a = 0.5 / gaP
                D[ii, jj, 0, 0, 1] = PB
                D[ii, jj, 1, 0, 1] = a
                D[ii, jj, 0, 1, 0] = PA
                D[ii, jj, 1, 2, 0] = a
                D[ii, jj, 0, 1, 1] = PB * D[ii, jj, 0, 1, 0] + D[ii, jj, 1, 2, 0]
                D[ii, jj, 1, 2, 1] = a * D[ii, jj, 0, 1, 0] + PB * D[ii, jj, 1, 2, 0]
                D[ii, jj, 2, 1, 1] = a * D[ii, jj, 1, 2, 0]

                D[ii, jj, 0, 2, 0] = PA * D[ii, jj, 0, 1, 0] + D[ii, jj, 1, 2, 0]
                D[ii, jj, 1, 2, 0] = a * D[ii, jj, 0, 1, 0] + PA * D[ii, jj, 1, 2, 0]
                D[ii, jj, 2, 2, 0] = a * D[ii, jj, 1, 2, 0]
                D[ii, jj, 0, 2, 1] = PB * D[ii, jj, 0, 2, 0] + D[ii, jj, 1, 2, 0]
                D[ii, jj, 1, 2, 1] = a * D[ii, jj, 0, 2, 0] + PB * D[ii, jj, 1, 2, 0] + 2. * D[ii, jj, 2, 2, 0]
                D[ii, jj, 2, 2, 1] = a * D[ii, jj, 1, 2, 0] + PB * D[ii, jj, 2, 2, 0]
                D[ii, jj, 3, 2, 1] = a * D[ii, jj, 2, 2, 0]
                D[ii, jj, 0, 3, 0] = PA * D[ii, jj, 0, 2, 0] + D[ii, jj, 1, 2, 0]
                D[ii, jj, 1, 3, 0] = a * D[ii, jj, 0, 2, 0] + PA * D[ii, jj, 1, 2, 0] + 2. * D[ii, jj, 2, 2, 0]
                D[ii, jj, 2, 3, 0] = a * D[ii, jj, 1, 2, 0] + PA * D[ii, jj, 2, 2, 0]
                D[ii, jj, 3, 3, 0] = a * D[ii, jj, 2, 2, 0]
                D[ii, jj, 0, 3, 1] = PB * D[ii, jj, 0, 3, 0] + D[ii, jj, 1, 3, 0]
                D[ii, jj, 1, 3, 1] = a * D[ii, jj, 0, 3, 0] + PB * D[ii, jj, 1, 3, 0] + 2. * D[ii, jj, 2, 3, 0]
                D[ii, jj, 2, 3, 1] = a * D[ii, jj, 1, 3, 0] + PB * D[ii, jj, 2, 3, 0] + 3. * D[ii, jj, 3, 3, 0]
                D[ii, jj, 3, 3, 1] = a * D[ii, jj, 2, 3, 0] + PB * D[ii, jj, 3, 3, 0]
                D[ii, jj, 4, 3, 1] = a * D[ii, jj, 3, 3, 0]

                D[jj, ii, 0, 1, 3] = D[ii, jj, 0, 3, 1]
                D[jj, ii, 1, 2, 3] = D[ii, jj, 1, 3, 1]
                D[jj, ii, 2, 1, 3] = D[ii, jj, 2, 3, 1]
                D[jj, ii, 3, 1, 3] = D[ii, jj, 3, 3, 1]
                D[jj, ii, 4, 1, 3] = D[ii, jj, 4, 3, 1]
            elif l1 == 3 and l2 == 2:
                a = 0.5 / gaP
                D[ii, jj, 0, 0, 1] = PB
                D[ii, jj, 1, 0, 1] = a
                D[ii, jj, 0, 0, 2] = PB * D[ii, jj, 0, 0, 1] + D[ii, jj, 1, 0, 1]
                D[ii, jj, 1, 0, 2] = a * D[ii, jj, 0, 0, 1] + PB * D[ii, jj, 1, 0, 1]
                D[ii, jj, 2, 0, 2] = a * D[ii, jj, 1, 0, 1]
                D[ii, jj, 0, 1, 0] = PA
                D[ii, jj, 1, 2, 0] = a
                D[ii, jj, 0, 1, 1] = PB * D[ii, jj, 0, 1, 0] + D[ii, jj, 1, 2, 0]
                D[ii, jj, 1, 2, 1] = a * D[ii, jj, 0, 1, 0] + PB * D[ii, jj, 1, 2, 0]
                D[ii, jj, 2, 1, 1] = a * D[ii, jj, 1, 2, 0]
                D[ii, jj, 0, 1, 2] = PB * D[ii, jj, 0, 1, 1] + D[ii, jj, 1, 2, 1]
                D[ii, jj, 1, 2, 2] = a * D[ii, jj, 0, 1, 1] + PB * D[ii, jj, 1, 2, 1] + 2. * D[ii, jj, 2, 1, 1]
                D[ii, jj, 2, 1, 2] = a * D[ii, jj, 1, 2, 1] + PB * D[ii, jj, 2, 1, 1]
                D[ii, jj, 3, 1, 2] = a * D[ii, jj, 2, 1, 1]
                D[ii, jj, 0, 2, 0] = PA * D[ii, jj, 0, 1, 0] + D[ii, jj, 1, 2, 0]
                D[ii, jj, 1, 2, 0] = a * D[ii, jj, 0, 1, 0] + PA * D[ii, jj, 1, 2, 0]
                D[ii, jj, 2, 2, 0] = a * D[ii, jj, 1, 2, 0]
                D[ii, jj, 0, 2, 1] = PB * D[ii, jj, 0, 2, 0] + D[ii, jj, 1, 2, 0]
                D[ii, jj, 1, 2, 1] = a * D[ii, jj, 0, 2, 0] + PB * D[ii, jj, 1, 2, 0] + 2. * D[ii, jj, 2, 2, 0]
                D[ii, jj, 2, 2, 1] = a * D[ii, jj, 1, 2, 0] + PB * D[ii, jj, 2, 2, 0]
                D[ii, jj, 3, 2, 1] = a * D[ii, jj, 2, 2, 0]
                D[ii, jj, 0, 2, 2] = PB * D[ii, jj, 0, 2, 1] + D[ii, jj, 1, 2, 1]
                D[ii, jj, 1, 2, 2] = a * D[ii, jj, 0, 2, 1] + PB * D[ii, jj, 1, 2, 1] + 2. * D[ii, jj, 2, 2, 1]
                D[ii, jj, 2, 2, 2] = a * D[ii, jj, 1, 2, 1] + PB * D[ii, jj, 2, 2, 1] + 3. * D[ii, jj, 3, 2, 1]
                D[ii, jj, 3, 2, 2] = a * D[ii, jj, 2, 2, 1] + PB * D[ii, jj, 3, 2, 1]
                D[ii, jj, 4, 2, 2] = a * D[ii, jj, 3, 2, 1]
                D[ii, jj, 0, 3, 0] = PA * D[ii, jj, 0, 2, 0] + D[ii, jj, 1, 2, 0]
                D[ii, jj, 1, 3, 0] = a * D[ii, jj, 0, 2, 0] + PA * D[ii, jj, 1, 2, 0] + 2. * D[ii, jj, 2, 2, 0]
                D[ii, jj, 2, 3, 0] = a * D[ii, jj, 1, 2, 0] + PA * D[ii, jj, 2, 2, 0]
                D[ii, jj, 3, 3, 0] = a * D[ii, jj, 2, 2, 0]
                D[ii, jj, 0, 3, 1] = PB * D[ii, jj, 0, 3, 0] + D[ii, jj, 1, 3, 0]
                D[ii, jj, 1, 3, 1] = a * D[ii, jj, 0, 3, 0] + PB * D[ii, jj, 1, 3, 0] + 2. * D[ii, jj, 2, 3, 0]
                D[ii, jj, 2, 3, 1] = a * D[ii, jj, 1, 3, 0] + PB * D[ii, jj, 2, 3, 0] + 3. * D[ii, jj, 3, 3, 0]
                D[ii, jj, 3, 3, 1] = a * D[ii, jj, 2, 3, 0] + PB * D[ii, jj, 3, 3, 0]
                D[ii, jj, 4, 3, 1] = a * D[ii, jj, 3, 3, 0]
                D[ii, jj, 0, 3, 2] = PB * D[ii, jj, 0, 3, 1] + D[ii, jj, 1, 3, 1]
                D[ii, jj, 1, 3, 2] = a * D[ii, jj, 0, 3, 1] + PB * D[ii, jj, 1, 3, 1] + 2. * D[ii, jj, 2, 3, 1]
                D[ii, jj, 2, 3, 2] = a * D[ii, jj, 1, 3, 1] + PB * D[ii, jj, 2, 3, 1] + 3. * D[ii, jj, 3, 3, 1]
                D[ii, jj, 3, 3, 2] = a * D[ii, jj, 2, 3, 1] + PB * D[ii, jj, 3, 3, 1] + 4. * D[ii, jj, 4, 3, 1]
                D[ii, jj, 4, 3, 2] = a * D[ii, jj, 3, 3, 1] + PB * D[ii, jj, 4, 3, 1]
                D[ii, jj, 5, 3, 2] = a * D[ii, jj, 4, 3, 1]

                D[jj, ii, 0, 2, 3] = D[ii, jj, 0, 3, 2]
                D[jj, ii, 1, 2, 3] = D[ii, jj, 1, 3, 2]
                D[jj, ii, 2, 2, 3] = D[ii, jj, 2, 3, 2]
                D[jj, ii, 3, 2, 3] = D[ii, jj, 3, 3, 2]
                D[jj, ii, 4, 2, 3] = D[ii, jj, 4, 3, 2]
                D[jj, ii, 5, 2, 3] = D[ii, jj, 5, 3, 2]
            elif l1 == 3 and l2 == 3:
                a = 0.5 / gaP
                D[ii, jj, 0, 0, 1] = PB
                D[ii, jj, 0, 0, 1] = PB
                D[ii, jj, 1, 0, 1] = a
                D[ii, jj, 0, 0, 2] = PB * D[ii, jj, 0, 0, 1] + D[ii, jj, 1, 0, 1]
                D[ii, jj, 1, 0, 2] = a * D[ii, jj, 0, 0, 1] + PB * D[ii, jj, 1, 0, 1]
                D[ii, jj, 2, 0, 2] = a * D[ii, jj, 1, 0, 1]
                D[ii, jj, 0, 0, 3] = PB * D[ii, jj, 0, 0, 2] + D[ii, jj, 1, 0, 2]
                D[ii, jj, 1, 0, 3] = a * D[ii, jj, 0, 0, 2] + PB * D[ii, jj, 1, 0, 2] + 2. * D[ii, jj, 2, 0, 2]
                D[ii, jj, 2, 0, 3] = a * D[ii, jj, 1, 0, 2] + PB * D[ii, jj, 2, 0, 2]
                D[ii, jj, 3, 0, 3] = a * D[ii, jj, 2, 0, 2]
                D[ii, jj, 0, 1, 0] = PA
                D[ii, jj, 1, 2, 0] = a
                D[ii, jj, 0, 1, 1] = PB * D[ii, jj, 0, 1, 0] + D[ii, jj, 1, 2, 0]
                D[ii, jj, 1, 2, 1] = a * D[ii, jj, 0, 1, 0] + PB * D[ii, jj, 1, 2, 0]
                D[ii, jj, 2, 1, 1] = a * D[ii, jj, 1, 2, 0]
                D[ii, jj, 0, 1, 2] = PB * D[ii, jj, 0, 1, 1] + D[ii, jj, 1, 2, 1]
                D[ii, jj, 1, 2, 2] = a * D[ii, jj, 0, 1, 1] + PB * D[ii, jj, 1, 2, 1] + 2. * D[ii, jj, 2, 1, 1]
                D[ii, jj, 2, 1, 2] = a * D[ii, jj, 1, 2, 1] + PB * D[ii, jj, 2, 1, 1]
                D[ii, jj, 3, 1, 2] = a * D[ii, jj, 2, 1, 1]

                D[ii, jj, 0, 1, 3] = PB * D[ii, jj, 0, 1, 2] + D[ii, jj, 1, 2, 2]
                D[ii, jj, 1, 2, 3] = a * D[ii, jj, 0, 1, 2] + PB * D[ii, jj, 1, 2, 2] + 2. * D[ii, jj, 2, 1, 2]
                D[ii, jj, 2, 1, 3] = a * D[ii, jj, 1, 2, 2] + PB * D[ii, jj, 2, 1, 2] + 3. * D[ii, jj, 3, 1, 2]
                D[ii, jj, 3, 1, 3] = a * D[ii, jj, 2, 1, 2] + PB * D[ii, jj, 3, 1, 2]
                D[ii, jj, 4, 1, 3] = a * D[ii, jj, 3, 1, 2]
                D[ii, jj, 0, 2, 0] = PA * D[ii, jj, 0, 1, 0] + D[ii, jj, 1, 2, 0]
                D[ii, jj, 1, 2, 0] = a * D[ii, jj, 0, 1, 0] + PA * D[ii, jj, 1, 2, 0]
                D[ii, jj, 2, 2, 0] = a * D[ii, jj, 1, 2, 0]
                D[ii, jj, 0, 2, 1] = PB * D[ii, jj, 0, 2, 0] + D[ii, jj, 1, 2, 0]
                D[ii, jj, 1, 2, 1] = a * D[ii, jj, 0, 2, 0] + PB * D[ii, jj, 1, 2, 0] + 2. * D[ii, jj, 2, 2, 0]
                D[ii, jj, 2, 2, 1] = a * D[ii, jj, 1, 2, 0] + PB * D[ii, jj, 2, 2, 0]
                D[ii, jj, 3, 2, 1] = a * D[ii, jj, 2, 2, 0]
                D[ii, jj, 0, 2, 2] = PB * D[ii, jj, 0, 2, 1] + D[ii, jj, 1, 2, 1]
                D[ii, jj, 1, 2, 2] = a * D[ii, jj, 0, 2, 1] + PB * D[ii, jj, 1, 2, 1] + 2. * D[ii, jj, 2, 2, 1]
                D[ii, jj, 2, 2, 2] = a * D[ii, jj, 1, 2, 1] + PB * D[ii, jj, 2, 2, 1] + 3. * D[ii, jj, 3, 2, 1]
                D[ii, jj, 3, 2, 2] = a * D[ii, jj, 2, 2, 1] + PB * D[ii, jj, 3, 2, 1]
                D[ii, jj, 4, 2, 2] = a * D[ii, jj, 3, 2, 1]

                D[ii, jj, 0, 2, 3] = PB * D[ii, jj, 0, 2, 2] + D[ii, jj, 1, 2, 2]
                D[ii, jj, 1, 2, 3] = a * D[ii, jj, 0, 2, 2] + PB * D[ii, jj, 1, 2, 2] + 2. * D[ii, jj, 2, 2, 2]
                D[ii, jj, 2, 2, 3] = a * D[ii, jj, 1, 2, 2] + PB * D[ii, jj, 2, 2, 2] + 3. * D[ii, jj, 3, 2, 2]
                D[ii, jj, 3, 2, 3] = a * D[ii, jj, 2, 2, 2] + PB * D[ii, jj, 3, 2, 2] + 4. * D[ii, jj, 4, 2, 2]
                D[ii, jj, 4, 2, 3] = a * D[ii, jj, 3, 2, 2] + PB * D[ii, jj, 4, 2, 2]
                D[ii, jj, 5, 2, 3] = a * D[ii, jj, 4, 2, 2]
                D[ii, jj, 0, 3, 0] = PA * D[ii, jj, 0, 2, 0] + D[ii, jj, 1, 2, 0]
                D[ii, jj, 1, 3, 0] = a * D[ii, jj, 0, 2, 0] + PA * D[ii, jj, 1, 2, 0] + 2. * D[ii, jj, 2, 2, 0]
                D[ii, jj, 2, 3, 0] = a * D[ii, jj, 1, 2, 0] + PA * D[ii, jj, 2, 2, 0]
                D[ii, jj, 3, 3, 0] = a * D[ii, jj, 2, 2, 0]
                D[ii, jj, 0, 3, 1] = PB * D[ii, jj, 0, 3, 0] + D[ii, jj, 1, 3, 0]
                D[ii, jj, 1, 3, 1] = a * D[ii, jj, 0, 3, 0] + PB * D[ii, jj, 1, 3, 0] + 2. * D[ii, jj, 2, 3, 0]
                D[ii, jj, 2, 3, 1] = a * D[ii, jj, 1, 3, 0] + PB * D[ii, jj, 2, 3, 0] + 3. * D[ii, jj, 3, 3, 0]
                D[ii, jj, 3, 3, 1] = a * D[ii, jj, 2, 3, 0] + PB * D[ii, jj, 3, 3, 0]
                D[ii, jj, 4, 3, 1] = a * D[ii, jj, 3, 3, 0]

                D[ii, jj, 0, 3, 2] = PB * D[ii, jj, 0, 3, 1] + D[ii, jj, 1, 3, 1]
                D[ii, jj, 1, 3, 2] = a * D[ii, jj, 0, 3, 1] + PB * D[ii, jj, 1, 3, 1] + 2. * D[ii, jj, 2, 3, 1]
                D[ii, jj, 2, 3, 2] = a * D[ii, jj, 1, 3, 1] + PB * D[ii, jj, 2, 3, 1] + 3. * D[ii, jj, 3, 3, 1]
                D[ii, jj, 3, 3, 2] = a * D[ii, jj, 2, 3, 1] + PB * D[ii, jj, 3, 3, 1] + 4. * D[ii, jj, 4, 3, 1]
                D[ii, jj, 4, 3, 2] = a * D[ii, jj, 3, 3, 1] + PB * D[ii, jj, 4, 3, 1]
                D[ii, jj, 5, 3, 2] = a * D[ii, jj, 4, 3, 1]

                D[ii, jj, 0, 3, 3] = PB * D[ii, jj, 0, 3, 2] + D[ii, jj, 1, 3, 2]
                D[ii, jj, 1, 3, 3] = a * D[ii, jj, 0, 3, 2] + PB * D[ii, jj, 1, 3, 2] + 2. * D[ii, jj, 2, 3, 2]
                D[ii, jj, 2, 3, 3] = a * D[ii, jj, 1, 3, 2] + PB * D[ii, jj, 2, 3, 2] + 3. * D[ii, jj, 3, 3, 2]
                D[ii, jj, 3, 3, 3] = a * D[ii, jj, 2, 3, 2] + PB * D[ii, jj, 3, 3, 2] + 4. * D[ii, jj, 4, 3, 2]
                D[ii, jj, 4, 3, 3] = a * D[ii, jj, 3, 3, 2] + PB * D[ii, jj, 4, 3, 2] + 5. * D[ii, jj, 5, 3, 2]
                D[ii, jj, 5, 3, 3] = a * D[ii, jj, 4, 3, 2] + PB * D[ii, jj, 5, 3, 2]
                D[ii, jj, 6, 3, 3] = a * D[ii, jj, 5, 3, 2]

                D[jj, ii, 0, 3, 3] = D[ii, jj, 0, 3, 3]
                D[jj, ii, 1, 3, 3] = D[ii, jj, 1, 3, 3]
                D[jj, ii, 2, 3, 3] = D[ii, jj, 2, 3, 3]
                D[jj, ii, 3, 3, 3] = D[ii, jj, 3, 3, 3]
                D[jj, ii, 4, 3, 3] = D[ii, jj, 4, 3, 3]
                D[jj, ii, 5, 3, 3] = D[ii, jj, 5, 3, 3]
                D[jj, ii, 6, 3, 3] = D[ii, jj, 6, 3, 3]

            else:
                print('case not programmed yet: l1/2= ', str(l1), str(l2))

    for i in range(0, N):
        j = i

        gaP = ga[i] + ga[j]
        Px = (ga[i] * x[i] + ga[j] * x[j]) / gaP
        PA = Px - x[i]
        PB = Px - x[j]

        l1 = l[i]
        l2 = l[j]

        if l1 == 0 and l2 == 0:
            D[i, j, 0, 0, 0] = 1
        elif l1 == 1:
            a = 0.5 / gaP
            D[i, j, 0, 0, 1] = PB
            D[i, j, 1, 0, 1] = a
            D[i, j, 0, 1, 0] = PA
            D[i, j, 1, 1, 0] = a  # corrected

            D[i, j, 0, 1, 1] = PB * D[i, j, 0, 1, 0] + D[i, j, 1, 1, 0]  # corrected
            D[i, j, 1, 1, 1] = a * D[i, j, 0, 1, 0] + PB * D[i, j, 1, 1, 0]  # corrected
            D[i, j, 2, 1, 1] = a * D[i, j, 1, 1, 0]  # corrected
        elif l1 == 2:
            a = 0.5 / gaP

            D[i, j, 0, 0, 1] = PB
            D[i, j, 1, 0, 1] = a
            D[i, j, 0, 0, 2] = PB * D[i, j, 0, 0, 1] + D[i, j, 1, 0, 1]
            D[i, j, 1, 0, 2] = a * D[i, j, 0, 0, 1] + PB * D[i, j, 1, 0, 1]
            D[i, j, 2, 0, 2] = a * D[i, j, 1, 0, 1]
            D[i, j, 0, 1, 0] = PA
            D[i, j, 1, 1, 0] = a  # corrected
            D[i, j, 0, 1, 1] = PB * D[i, j, 0, 1, 0] + D[i, j, 1, 1, 0]  # corrected
            D[i, j, 1, 1, 1] = a * D[i, j, 0, 1, 0] + PB * D[i, j, 1, 1, 0]  # corrected
            D[i, j, 2, 1, 1] = a * D[i, j, 1, 1, 0]  # corrected
            D[i, j, 0, 1, 2] = PB * D[i, j, 0, 1, 1] + D[i, j, 1, 1, 1]  # corrected
            D[i, j, 1, 1, 2] = a * D[i, j, 0, 1, 1] + PB * D[i, j, 1, 1, 1] + 2. * D[i, j, 2, 1, 1]  # corrected
            D[i, j, 2, 1, 2] = a * D[i, j, 1, 1, 1] + PB * D[i, j, 2, 1, 1]
            D[i, j, 3, 1, 2] = a * D[i, j, 2, 1, 1]

            D[i, j, 0, 2, 0] = PA * D[i, j, 0, 1, 0] + D[i, j, 1, 1, 0]  # corrected
            D[i, j, 1, 2, 0] = a * D[i, j, 0, 1, 0] + PA * D[i, j, 1, 1, 0]  # corrected
            D[i, j, 2, 2, 0] = a * D[i, j, 1, 1, 0]  # corrected
            D[i, j, 0, 2, 1] = PB * D[i, j, 0, 2, 0] + D[i, j, 1, 2, 0]
            D[i, j, 1, 2, 1] = a * D[i, j, 0, 2, 0] + PB * D[i, j, 1, 2, 0] + 2. * D[i, j, 2, 2, 0]
            D[i, j, 2, 2, 1] = a * D[i, j, 1, 2, 0] + PB * D[i, j, 2, 2, 0]
            D[i, j, 3, 2, 1] = a * D[i, j, 2, 2, 0]
            D[i, j, 0, 2, 2] = PB * D[i, j, 0, 2, 1] + D[i, j, 1, 2, 1]
            D[i, j, 1, 2, 2] = a * D[i, j, 0, 2, 1] + PB * D[i, j, 1, 2, 1] + 2. * D[i, j, 2, 2, 1]
            D[i, j, 2, 2, 2] = a * D[i, j, 1, 2, 1] + PB * D[i, j, 2, 2, 1] + 3. * D[i, j, 3, 2, 1]
            D[i, j, 3, 2, 2] = a * D[i, j, 2, 2, 1] + PB * D[i, j, 3, 2, 1]
            D[i, j, 4, 2, 2] = a * D[i, j, 3, 2, 1]
        elif l1 == 3:
            a = 0.5 / gaP
            D[i, j, 0, 0, 1] = PB
            D[i, j, 0, 0, 1] = PB
            D[i, j, 1, 0, 1] = a
            D[i, j, 0, 0, 2] = PB * D[i, j, 0, 0, 1] + D[i, j, 1, 0, 1]
            D[i, j, 1, 0, 2] = a * D[i, j, 0, 0, 1] + PB * D[i, j, 1, 0, 1]
            D[i, j, 2, 0, 2] = a * D[i, j, 1, 0, 1]
            D[i, j, 0, 0, 3] = PB * D[i, j, 0, 0, 2] + D[i, j, 1, 0, 2]
            D[i, j, 1, 0, 3] = a * D[i, j, 0, 0, 2] + PB * D[i, j, 1, 0, 2] + 2. * D[i, j, 2, 0, 2]
            D[i, j, 2, 0, 3] = a * D[i, j, 1, 0, 2] + PB * D[i, j, 2, 0, 2]
            D[i, j, 3, 0, 3] = a * D[i, j, 2, 0, 2]
            D[i, j, 0, 1, 0] = PA
            D[i, j, 1, 1, 0] = a  # corrected
            D[i, j, 0, 1, 1] = PB * D[i, j, 0, 1, 0] + D[i, j, 1, 1, 0]  # corrected
            D[i, j, 1, 1, 1] = a * D[i, j, 0, 1, 0] + PB * D[i, j, 1, 1, 0]  # corrected
            D[i, j, 2, 1, 1] = a * D[i, j, 1, 1, 0]  # corrected
            D[i, j, 0, 1, 2] = PB * D[i, j, 0, 1, 1] + D[i, j, 1, 1, 1]  # corrected
            D[i, j, 1, 1, 2] = a * D[i, j, 0, 1, 1] + PB * D[i, j, 1, 1, 1] + 2. * D[i, j, 2, 1, 1]  # corrected
            D[i, j, 2, 1, 2] = a * D[i, j, 1, 1, 1] + PB * D[i, j, 2, 1, 1]  # corrected
            D[i, j, 3, 1, 2] = a * D[i, j, 2, 1, 1]

            D[i, j, 0, 1, 3] = PB * D[i, j, 0, 1, 2] + D[i, j, 1, 1, 2]  # corrected
            D[i, j, 1, 1, 3] = a * D[i, j, 0, 1, 2] + PB * D[i, j, 1, 1, 2] + 2. * D[i, j, 2, 1, 2]
            D[i, j, 2, 1, 3] = a * D[i, j, 1, 1, 2] + PB * D[i, j, 2, 1, 2] + 3. * D[i, j, 3, 1, 2]
            D[i, j, 3, 1, 3] = a * D[i, j, 2, 1, 2] + PB * D[i, j, 3, 1, 2]
            D[i, j, 4, 1, 3] = a * D[i, j, 3, 1, 2]
            D[i, j, 0, 2, 0] = PA * D[i, j, 0, 1, 0] + D[i, j, 1, 1, 0]  # corrected
            D[i, j, 1, 2, 0] = a * D[i, j, 0, 1, 0] + PA * D[i, j, 1, 1, 0]  # corrected
            D[i, j, 2, 2, 0] = a * D[i, j, 1, 1, 0]  # corrected
            D[i, j, 0, 2, 1] = PB * D[i, j, 0, 2, 0] + D[i, j, 1, 2, 0]
            D[i, j, 1, 2, 1] = a * D[i, j, 0, 2, 0] + PB * D[i, j, 1, 2, 0] + 2. * D[i, j, 2, 2, 0]
            D[i, j, 2, 2, 1] = a * D[i, j, 1, 2, 0] + PB * D[i, j, 2, 2, 0]
            D[i, j, 3, 2, 1] = a * D[i, j, 2, 2, 0]
            D[i, j, 0, 2, 2] = PB * D[i, j, 0, 2, 1] + D[i, j, 1, 2, 1]
            D[i, j, 1, 2, 2] = a * D[i, j, 0, 2, 1] + PB * D[i, j, 1, 2, 1] + 2. * D[i, j, 2, 2, 1]
            D[i, j, 2, 2, 2] = a * D[i, j, 1, 2, 1] + PB * D[i, j, 2, 2, 1] + 3. * D[i, j, 3, 2, 1]
            D[i, j, 3, 2, 2] = a * D[i, j, 2, 2, 1] + PB * D[i, j, 3, 2, 1]
            D[i, j, 4, 2, 2] = a * D[i, j, 3, 2, 1]

            D[i, j, 0, 2, 3] = PB * D[i, j, 0, 2, 2] + D[i, j, 1, 2, 2]
            D[i, j, 1, 2, 3] = a * D[i, j, 0, 2, 2] + PB * D[i, j, 1, 2, 2] + 2. * D[i, j, 2, 2, 2]
            D[i, j, 2, 2, 3] = a * D[i, j, 1, 2, 2] + PB * D[i, j, 2, 2, 2] + 3. * D[i, j, 3, 2, 2]
            D[i, j, 3, 2, 3] = a * D[i, j, 2, 2, 2] + PB * D[i, j, 3, 2, 2] + 4. * D[i, j, 4, 2, 2]
            D[i, j, 4, 2, 3] = a * D[i, j, 3, 2, 2] + PB * D[i, j, 4, 2, 2]
            D[i, j, 5, 2, 3] = a * D[i, j, 4, 2, 2]
            D[i, j, 0, 3, 0] = PA * D[i, j, 0, 2, 0] + D[i, j, 1, 2, 0]
            D[i, j, 1, 3, 0] = a * D[i, j, 0, 2, 0] + PA * D[i, j, 1, 2, 0] + 2. * D[i, j, 2, 2, 0]
            D[i, j, 2, 3, 0] = a * D[i, j, 1, 2, 0] + PA * D[i, j, 2, 2, 0]
            D[i, j, 3, 3, 0] = a * D[i, j, 2, 2, 0]
            D[i, j, 0, 3, 1] = PB * D[i, j, 0, 3, 0] + D[i, j, 1, 3, 0]
            D[i, j, 1, 3, 1] = a * D[i, j, 0, 3, 0] + PB * D[i, j, 1, 3, 0] + 2. * D[i, j, 2, 3, 0]
            D[i, j, 2, 3, 1] = a * D[i, j, 1, 3, 0] + PB * D[i, j, 2, 3, 0] + 3. * D[i, j, 3, 3, 0]
            D[i, j, 3, 3, 1] = a * D[i, j, 2, 3, 0] + PB * D[i, j, 3, 3, 0]
            D[i, j, 4, 3, 1] = a * D[i, j, 3, 3, 0]

            D[i, j, 0, 3, 2] = PB * D[i, j, 0, 3, 1] + D[i, j, 1, 3, 1]
            D[i, j, 1, 3, 2] = a * D[i, j, 0, 3, 1] + PB * D[i, j, 1, 3, 1] + 2. * D[i, j, 2, 3, 1]
            D[i, j, 2, 3, 2] = a * D[i, j, 1, 3, 1] + PB * D[i, j, 2, 3, 1] + 3. * D[i, j, 3, 3, 1]
            D[i, j, 3, 3, 2] = a * D[i, j, 2, 3, 1] + PB * D[i, j, 3, 3, 1] + 4. * D[i, j, 4, 3, 1]
            D[i, j, 4, 3, 2] = a * D[i, j, 3, 3, 1] + PB * D[i, j, 4, 3, 1]
            D[i, j, 5, 3, 2] = a * D[i, j, 4, 3, 1]

            D[i, j, 0, 3, 3] = PB * D[i, j, 0, 3, 2] + D[i, j, 1, 3, 2]
            D[i, j, 1, 3, 3] = a * D[i, j, 0, 3, 2] + PB * D[i, j, 1, 3, 2] + 2. * D[i, j, 2, 3, 2]
            D[i, j, 2, 3, 3] = a * D[i, j, 1, 3, 2] + PB * D[i, j, 2, 3, 2] + 3. * D[i, j, 3, 3, 2]
            D[i, j, 3, 3, 3] = a * D[i, j, 2, 3, 2] + PB * D[i, j, 3, 3, 2] + 4. * D[i, j, 4, 3, 2]
            D[i, j, 4, 3, 3] = a * D[i, j, 3, 3, 2] + PB * D[i, j, 4, 3, 2] + 5. * D[i, j, 5, 3, 2]
            D[i, j, 5, 3, 3] = a * D[i, j, 4, 3, 2] + PB * D[i, j, 5, 3, 2]
            D[i, j, 6, 3, 3] = a * D[i, j, 5, 3, 2]

        else:
            print('case not programmed yet: l1/2= ', str(l1), str(l2))
    return D
