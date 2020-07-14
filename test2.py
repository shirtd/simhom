from __future__ import annotations

from simhom.sequences import *

if __name__ == '__main__':

    # KS = [[3, 0, 1], [2, 0, 1]]
    # LS = [[3, 0, 1], [2, 0, 1]]

    KS = [[0,1,3],
        [1,3,4],
        [1,2,4],
        [2,4,5],
        [0,2,5],
        [0,3,5]]

    LS = [[0,1],[1,2],[0,2],
        [3,4],[4,5],[3,5]]

    K = SimplicialComplex(list(map(Simplex, KS)), 'K_1')
    L = SimplicialComplex(list(map(Simplex, LS)), 'L_1')

    S = ChainComplexPairSES.simplicial_init(K, L)

    H = HomologyPairLES.complex_init(S.cycles() / S.boundaries())

    # # S = [[0,1,2,3], [3,4],[2,4], [5,6,7], [6,7,8], [5,7,8], [5,6,8]]
    # # # S = [[0,1,2], [3,4],[4,5],[3,5]]
    # # # # S = [[0,1],[1,2],[0,2],[1,2,3]]
    # # # # S = [[0,1], [2,3]]
    # # K = SimplicialComplex(lmap(Simplex, S))
    #
    # # KS = [
    # #     [0,1,2,3],
    # #     [1,2,3,4],
    # #     [0,1,3,5],
    # #     [1,3,4,5]
    # #     # [0,3,5],
    # #     # [3,4,5],
    # #     # [1,4,5],
    # #     # [0,1,5]
    # # ] # [[0,1,2], [1,2,3]]
    # # LS = [
    # #     [0,3,5],
    # #     [3,4,5],
    # #     [1,4,5],
    # #     [0,1,5],
    # #     [0,2,1],
    # #     [0,3,2],
    # #     [1,2,4],
    # #     [2,3,4]
    # # ] # [[0,1],[0,2],[2,3],[1,3]]
    #
    # KS = [[0,1,2],[0,1,3],[1,2,4],
    #     [0,2],[0,5],[2,5],
    #     [0,1,5],[1,2,5]]
    # LS = [[0,1,3],[1,2,4],
    #     [0,2],[0,5],[2,5]]
    #
    # K = SimplicialComplex(list(map(Simplex, KS)), 'K_1')
    # L = SimplicialComplex(list(map(Simplex, LS)), 'L_1')
    #
    # C = SimplicialChainComplex.simplicial_init(K)
    # CL = SimplicialChainComplex.simplicial_init(L)
    #
    # S = ChainComplexPairSES.simplicial_init(K, L)
    # H = HomologyPairLES.complex_init(S.cycles() / S.boundaries())
    #
    # KS2 = [[0,1,2],[0,1,3],[1,2,4],
    #     [0,2],[0,5],[2,5],
    #     [0,1,5],[1,2,5]]
    # LS2 = [[0,1,3],[1,2,4],[0,2,5]]
    #
    # K2 = SimplicialComplex(list(map(Simplex, KS2)), 'K_2')
    # L2 = SimplicialComplex(list(map(Simplex, LS2)), 'L_2')
    #
    # C2 = SimplicialChainComplex.simplicial_init(K2)
    # CL2 = SimplicialChainComplex.simplicial_init(L2)
    #
    # S2 = ChainComplexPairSES.simplicial_init(K2, L2)
    # H2 = HomologyPairLES.complex_init(S2.cycles() / S2.boundaries())
