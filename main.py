from simhom.simplicial import SimplicialComplex
from simhom.homology import ChainGroup, Boundary
from simhom.util import lmap

if __name__ == '__main__':
    S = [[0,1,2],[3,4],[3,5],[5,6],[4,6],[7,8,9],[8,9,10],[7,9,10],[7,8,10]]
    K = SimplicialComplex(S)
    zero = ChainGroup(K, -1)
    C0, C1, C2 = lmap(lambda d: ChainGroup(K, d), (0,1,2))
    D0, D1, D2 = Boundary(C0, ChainGroup(K, -1)), Boundary(C1, C0), Boundary(C2, C1)

    H0 = D0.kernel() / D1.im
    H1 = D1.kernel() / D2.im
    H2 = D2.kernel() / ChainGroup(K, 2, [])
