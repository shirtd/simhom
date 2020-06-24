from simhom.simplicial import Simplex, SimplicialComplex, Chain, ZeroChain
from simhom.algebra import FGCoset, Quotient
from simhom.homology import ChainGroup, Boundary, ChainQuotient, ChainCoset, ChainComplex, Homology, ZeroChainCoset
# from simhom.homology import ChainGroup, Boundary, RelativeChainGroup, RelativeBoundary
from simhom.util import *



if __name__ == '__main__':
    # S = [[0,1,2,3], [3,4],[2,4], [5,6,7],[6,7,8],[5,7,8],[5,6,8]]
    # K = SimplicialComplex(lmap(Simplex, S))
    # C = ChainComplex(K)

    K = SimplicialComplex([Simplex(0,1,2),Simplex(0,2,3)])
    L = SimplicialComplex(lmap(Simplex, [[0,1],[1,2],[2,3],[0,3]]))
    CL = ChainComplex(L)
    CK = ChainComplex(K)
    C = ChainComplex(K, CL)

    H = Homology(C)
