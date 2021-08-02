import numpy as np


def pzerotablenew(L, M, N, nu):
    A = [L, M, N]
    A = np.sort(A)

    if L % 2 != 0 or M % 2 != 0 or N % 2 != 0:
        f2 = np.zeros((1, len(nu)))
    elif sum(A) == 0:
        f2 = np.ones((1, len(nu)))
    elif A[0] == 0 and A[1] == 0 and A[2] == 2:
        f2 = -(nu) ** 2 / 3
    elif A[0] == 0 and A[1] == 2 and A[2] == 2:
        f2 = ((nu) ** 4 / 15)
    elif A[0] == 2 and A[1] == 2 and A[2] == 2:
        f2 = -((nu) ** 6 / 105)
    elif A[0] == 0 and A[1] == 0 and A[2] == 4:
        f2 = ((nu) ** 4 / 5)
    elif A[0] == 0 and A[1] == 2 and A[2] == 4:
        f2 = -((nu) ** 6 / 35)
    elif A[0] == 2 and A[1] == 2 and A[2] == 4:
        f2 = ((nu) ** 8 / 315)
    elif A[0] == 2 and A[1] == 4 and A[2] == 4:
        f2 = -((nu) ** 10 / 1155)
    elif A[0] == 4 and A[1] == 4 and A[2] == 4:
        f2 = ((nu) ** 12 / 5005)
    elif A[0] == 0 and A[1] == 0 and A[2] == 6:
        f2 = -((nu) ** 6 / 7)
    elif A[0] == 0 and A[1] == 2 and A[2] == 6:
        f2 = ((nu) ** 8 / 63)
    elif A[0] == 2 and A[1] == 2 and A[2] == 6:
        f2 = -((nu) ** 10 / 693)
    elif A[0] == 2 and A[1] == 4 and A[2] == 6:
        f2 = ((nu) ** 12 / 3003)
    elif A[0] == 0 and A[1] == 4 and A[2] == 6:
        f2 = -((nu) ** 10 / 231)
    elif A[0] == 4 and A[1] == 4 and A[2] == 6:
        f2 = -((nu) ** 14 / 15015)
    elif A[0] == 2 and A[1] == 6 and A[2] == 6:
        f2 = -((nu) ** 14 / 9009)
    elif A[0] == 4 and A[1] == 6 and A[2] == 6:
        f2 = ((nu) ** 16 / 51051)
    elif A[0] == 6 and A[1] == 6 and A[2] == 6:
        f2 = -5 * ((nu) ** 18 / 969969)
    elif A[0] == 0 and A[1] == 6 and A[2] == 6:
        f2 = 5 * ((nu) ** 12 / 3003)
    elif A[0] == 0 and A[1] == 4 and A[2] == 4:
        f2 = ((nu) ** 8 / 105)
    elif A[0] == 0 and A[1] == 0 and A[2] == 8:
        f2 = ((nu) ** 8 / 9)
    elif A[0] == 0 and A[1] == 2 and A[2] == 8:
        f2 = -((nu) ** 10 / 99)
    elif A[0] == 2 and A[1] == 2 and A[2] == 8:
        f2 = ((nu) ** 12 / 1287)
    elif A[0] == 0 and A[1] == 4 and A[2] == 8:
        f2 = ((nu) ** 12 / 429)
    elif A[0] == 2 and A[1] == 4 and A[2] == 8:
        f2 = -((nu) ** 14 / 6435)
    elif A[0] == 4 and A[1] == 4 and A[2] == 8:
        f2 = ((nu) ** 16 / 36465)
    elif A[0] == 0 and A[1] == 6 and A[2] == 8:
        f2 = -((nu) ** 14 / 1287)
    elif A[0] == 2 and A[1] == 6 and A[2] == 8:
        f2 = ((nu) ** 16 / 21879)
    elif A[0] == 4 and A[1] == 6 and A[2] == 8:
        f2 = -((nu) ** 18 / 138567)
    elif A[0] == 6 and A[1] == 6 and A[2] == 8:
        f2 = 5 * ((nu) ** 20 / 2909907)
    elif A[0] == 0 and A[1] == 8 and A[2] == 8:
        f2 = 7 * ((nu) ** 16 / 21879)
    elif A[0] == 2 and A[1] == 8 and A[2] == 8:
        f2 = -7 * ((nu) ** 18 / 415701)
    elif A[0] == 4 and A[1] == 8 and A[2] == 8:
        f2 = ((nu) ** 20 / 415701)
    elif A[0] == 6 and A[1] == 8 and A[2] == 8:
        f2 = -5 * ((nu) ** 22 / 9561123)
    elif A[0] == 8 and A[1] == 8 and A[2] == 8:
        f2 = 7 * ((nu) ** 24 / 47805615)
    elif A[0] == 0 and A[1] == 0 and A[2] == 10:
        f2 = -((nu) ** 10 / 11)
    elif A[0] == 0 and A[1] == 2 and A[2] == 10:
        f2 = ((nu) ** 12 / 143)
    elif A[0] == 0 and A[1] == 4 and A[2] == 10:
        f2 = -((nu) ** 14 / 715)
    elif A[0] == 0 and A[1] == 6 and A[2] == 10:
        f2 = ((nu) ** 16 / 2431)
    elif A[0] == 0 and A[1] == 8 and A[2] == 10:
        f2 = -7 * ((nu) ** 18 / 46189)
    elif A[0] == 0 and A[1] == 10 and A[2] == 10:
        f2 = 3 * ((nu) ** 20 / 46189)
    elif A[0] == 0 and A[1] == 0 and A[2] == 12:
        f2 = ((nu) ** 12 / 13)
    elif A[0] == 0 and A[1] == 2 and A[2] == 12:
        f2 = -((nu) ** 14 / 195)
    elif A[0] == 0 and A[1] == 4 and A[2] == 12:
        f2 = ((nu) ** 16 / 1105)
    elif A[0] == 0 and A[1] == 6 and A[2] == 12:
        f2 = -((nu) ** 18 / 4199)
    elif A[0] == 0 and A[1] == 8 and A[2] == 12:
        f2 = ((nu) ** 20 / 12597)
    elif A[0] == 0 and A[1] == 10 and A[2] == 12:
        f2 = -3 * ((nu) ** 22 / 96577)
    elif A[0] == 0 and A[1] == 12 and A[2] == 12:
        f2 = 33 * ((nu) ** 24 / 2414425)
    elif A[0] == 0 and A[1] == 10 and A[2] == 10:
        f2 = 3 * ((nu) ** 20 / 46189)
    elif A[0] == 2 and A[1] == 10 and A[2] == 10:
        f2 = -3 * ((nu) ** 22 / 1062347)
    elif A[0] == 4 and A[1] == 10 and A[2] == 10:
        f2 = 9 * ((nu) ** 24 / 26558675)
    elif A[0] == 6 and A[1] == 10 and A[2] == 10:
        f2 = -((nu) ** 26 / 15935205)
    elif A[0] == 8 and A[1] == 10 and A[2] == 10:
        f2 = 7 * ((nu) ** 28 / 462120945)
    elif A[0] == 10 and A[1] == 10 and A[2] == 10:
        f2 = -21 * ((nu) ** 30 / 4775249765)
    elif A[0] == 2 and A[1] == 12 and A[2] == 12:
        f2 = -11 * ((nu) ** 26 / 21729825)
    elif A[0] == 4 and A[1] == 12 and A[2] == 12:
        f2 = 11 * ((nu) ** 28 / 210054975)
    elif A[0] == 6 and A[1] == 12 and A[2] == 12:
        f2 = -11 * ((nu) ** 30 / 1302340845)
    elif A[0] == 8 and A[1] == 12 and A[2] == 12:
        f2 = 7 * ((nu) ** 32 / 3907022535)
    elif A[0] == 10 and A[1] == 12 and A[2] == 12:
        f2 = -((nu) ** 34 / 2170568075)
    elif A[0] == 12 and A[1] == 12 and A[2] == 12:
        f2 = 11 * ((nu) ** 36 / 80311018775)
    elif A[0] == 2 and A[1] == 10 and A[2] == 12:
        f2 = 3 * ((nu) ** 24 / 2414425)
    elif A[0] == 4 and A[1] == 10 and A[2] == 12:
        f2 = -((nu) ** 26 / 7243275)
    elif A[0] == 6 and A[1] == 10 and A[2] == 12:
        f2 = (nu ** 28 / 42010995)
    elif A[0] == 8 and A[1] == 10 and A[2] == 12:
        f2 = -7 * ((nu) ** 30 / 1302340845)
    elif A[0] == 10 and A[1] == 10 and A[2] == 12:
        f2 = 7 * ((nu) ** 32 / 4775249765)
    elif A[0] == 2 and A[1] == 8 and A[2] == 12:
        f2 = -((nu) ** 22 / 289731)
    elif A[0] == 4 and A[1] == 8 and A[2] == 12:
        f2 = ((nu) ** 24 / 2414425)
    elif A[0] == 6 and A[1] == 8 and A[2] == 12:
        f2 = -((nu) ** 26 / 13037895)
    elif A[0] == 8 and A[1] == 8 and A[2] == 12:
        f2 = 7 * ((nu) ** 28 / 378098955)
    elif A[0] == 2 and A[1] == 6 and A[2] == 12:
        f2 = ((nu) ** 20 / 88179)
    elif A[0] == 4 and A[1] == 6 and A[2] == 12:
        f2 = -((nu) ** 22 / 676039)
    elif A[0] == 6 and A[1] == 6 and A[2] == 12:
        f2 = ((nu) ** 24 / 3380195)
    elif A[0] == 2 and A[1] == 4 and A[2] == 12:
        f2 = -((nu) ** 18 / 20995)
    elif A[0] == 4 and A[1] == 4 and A[2] == 12:
        f2 = ((nu) ** 20 / 146965)
    elif A[0] == 2 and A[1] == 2 and A[2] == 12:
        f2 = ((nu) ** 16 / 3315)
    elif A[0] == 2 and A[1] == 8 and A[2] == 10:
        f2 = ((nu) ** 20 / 138567)
    elif A[0] == 4 and A[1] == 8 and A[2] == 10:
        f2 = -((nu) ** 22 / 1062347)
    elif A[0] == 6 and A[1] == 8 and A[2] == 10:
        f2 = ((nu) ** 24 / 5311735)
    elif A[0] == 8 and A[1] == 8 and A[2] == 10:
        f2 = -7 * ((nu) ** 26 / 143416845)
    elif A[0] == 2 and A[1] == 6 and A[2] == 10:
        f2 = -((nu) ** 18 / 46189)
    elif A[0] == 4 and A[1] == 6 and A[2] == 10:
        f2 = ((nu) ** 20 / 323323)
    elif A[0] == 6 and A[1] == 6 and A[2] == 10:
        f2 = -5 * ((nu) ** 22 / 7436429)
    elif A[0] == 2 and A[1] == 4 and A[2] == 10:
        f2 = ((nu) ** 16 / 12155)
    elif A[0] == 4 and A[1] == 4 and A[2] == 10:
        f2 = -3 * (nu ** 18 / 230945)
    elif A[0] == 2 and A[1] == 2 and A[2] == 10:
        f2 = -(nu ** 14 / 2145)
    else:
        print('case not programmed yet')
        print([L, M, N])
        f2 = 0
    return f2
