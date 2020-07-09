from __future__ import annotations

from simhom.sequences import *

if __name__ == '__main__':
    # S = [[0,1,2,3], [3,4],[2,4], [5,6,7], [6,7,8], [5,7,8], [5,6,8]]
    # # S = [[0,1,2], [3,4],[4,5],[3,5]]
    # # # S = [[0,1],[1,2],[0,2],[1,2,3]]
    # # # S = [[0,1], [2,3]]
    # K = SimplicialComplex(lmap(Simplex, S))

    K = SimplicialComplex(list(map(Simplex, [[0,1,2], [1,2,3]])))
    L = SimplicialComplex(list(map(Simplex, [[0,1],[0,2],[2,3],[1,3]])))

    C = SimplicialChainComplex.simplicial_init(K)
    CL = SimplicialChainComplex.simplicial_init(L)

    # R = C / CL
    # H = R.cycles() / R.boundaries()

    S = ChainComplexPairSES.simplicial_init(K, L)

    H = S.cycles() / S.boundaries()
