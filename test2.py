from __future__ import annotations


from simhom.homology import *
from simhom.util import lmap, ker, im
import numpy as np

''' TODO ConnectingHomomorphism '''

if __name__ == '__main__':
    # # S = [[0,1,2,3], [3,4],[2,4], [5,6,7]]#, [5,6,7],[6,7,8],[5,7,8],[5,6,8]]
    # S = [[0,1,2], [3,4],[4,5],[3,5]]
    # # S = [[0,1],[1,2],[0,2],[1,2,3]]
    # # S = [[0,1], [2,3]]
    # K = SimplicialComplex(lmap(Simplex, S))
    K = SimplicialComplex(list(map(Simplex, [[0,1,2], [1,2,3]])))
    L = SimplicialComplex(list(map(Simplex, [[0,1],[0,2],[2,3],[1,3]])))


    CL = {d : ChainGroup.simplicial_init(L, d) for d in range(4)}

    HL2 = (CL[2] >> CL[1]).ker / (CL[3] >> CL[2]).im
    HL1 = (CL[1] >> CL[0]).ker / (CL[2] >> CL[1]).im
    HL0 = CL[0] / (CL[1] >> CL[0]).im

    print('HL2\t', HL2)
    print('HL1\t', HL1)
    print('HL0\t', HL0)


    C = {d : ChainGroup.simplicial_init(K, d) for d in range(4)}

    H2 = (C[2] >> C[1]).ker / (C[3] >> C[2]).im
    H1 = (C[1] >> C[0]).ker / (C[2] >> C[1]).im
    H0 = C[0] / (C[1] >> C[0]).im

    print('\nH2\t', H2)
    print('H1\t', H1)
    print('H0\t', H0)


    CR = {d : ChainQuotient(C[d], CL[d]) for d in range(4)}

    HR2 = (CR[2] >> CR[1]).ker / (CR[3] >> CR[2]).im
    HR1 = (CR[1] >> CR[0]).ker / (CR[2] >> CR[1]).im
    HR0 = CR[0] / (CR[1] >> CR[0]).im

    print('\nHR2\t', HR2)
    print('HR1\t', HR1)
    print('HR0\t', HR0)

    print(HL2 >> H2, H2 >> HR2)
    print(HL1 >> H1, H1 >> HR1)
    print(HL0 >> H0, H0 >> HR0)
