# symmetry.py

def findirrep(nirreps,slaterdet):

    dpta = [[1, 2], [2, 1]]
    dptb = [[1, 2, 3, 4], [2, 1, 4, 3], [3, 4, 1, 2], [4, 3, 2, 1]]
    dptc = [[1, 2, 3, 4, 5, 6, 7, 8], [2, 1, 4, 3, 6, 5, 8, 7], [3, 4, 1, 2, 7, 8, 5, 6],
            [4, 3, 2, 1, 8, 7, 6, 5], [5, 6, 7, 8, 1, 2, 3, 4], [6, 5, 8, 7, 2, 1, 4, 3],
            [7, 8, 5, 6, 3, 4, 1, 2], [8, 7, 6, 5, 4, 3, 2, 1]]
   
    if nirreps == 2:
        dpt = dpta
    elif nirreps == 4:
        dpt = dptb
    elif nirreps == 8:
        dpt = dptc

    sirr = []

    for irr in range(nirreps):
        occstr = slaterdet[irr].replace('0','').replace('a','1').replace('b','1')
        val = 0
        for digit in occstr:
            val = val + int(digit)
        if irr > 0 and val % 2 != 0:
            sirr.append(irr)

    nsirr = len(sirr)

    if nsirr == 0:
        irrep = 1
    elif nsirr == 1:
        irrep = sirr[0] + 1
    else:
        irrep = dpt[sirr[0]][sirr[1]]
        if nsirr > 2:
            for irr in range(2,nsirr):
                irrep = dpt[(irrep-1)][sirr[irr]]

    return irrep

