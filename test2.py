from __future__ import annotations


from simhom.chains import *
from simhom.util import hstack, lmap, cref
import numpy as np


def coset_reduce(a, H):
    if not len(H):
        return ((a, ZeroChain(a.dim)),)
    A = {a}
    for h in H:
        for x in A.copy():
            if any(s in h for s in x):
                A.add(min(x + h, x - h))
    return min(A)


class ZeroChainCoset(ZeroElement, AbstractChain):
    @classmethod
    def is_zero(cls, chain, subgroup):
        return chain in subgroup
    def __init__(self, chain, subgroup):
        ZeroElement.__init__(self)
        self._subgroup = subgroup
        self._rep = coset_reduce(chain, subgroup)
        AbstractChain.__init__(self, (), self._rep.dim)
    def __contains__(self, other):
        return other - self._rep in self._subgroup
    def __eq__(self, other):
        return ZeroElement.__eq__(self, other)
    def __hash__(self) -> int:
        return hash(self.id)
    def __mul__(self, a: int) -> AbstractChain:
        return self
    def __repr__(self) -> str:
        return '[%s[]]' % AbstractChain.__repr__(self)


class ChainCoset(GroupElement, AbstractChain):
    zero_t = ZeroChainCoset
    @classmethod
    def _mkzero(cls, chain, subgroup):
        return ZeroChainCoset(chain, subgroup)
    def __init__(self, chain, subgroup):
        self._subgroup = subgroup
        self._rep = coset_reduce(chain, subgroup)
        AbstractChain.__init__(self, self._rep.id, self._rep.dim)
    def __contains__(self, other):
        return other - self._rep in self._subgroup
    def __add__(self, other):
        return ChainCoset(self._rep + other._rep, self._subgroup)
    def __invert__(self):
        return self
    def __mul__(self, a):
        return ChainCoset(self._rep * a, self._subgroup)
    def __repr__(self):
        return '[%s]' % str(self._rep)
    def __lt__(self, other) -> bool:
        return self._rep < other._rep
    def __eq__(self, h) -> bool:
        return AbstractElement.__eq__(self, h)
    def __hash__(self) -> int:
        return AbstractElement.__hash__(self)


class ChainQuotient(Group[AbstractChain]):#, Basis[Chain]):
    def __init__(self, group, subgroup) -> None:
        self._group, self._subgroup, self._dim = group, subgroup, group.dim
        zero = ZeroChainCoset(group._zero, subgroup)
        Group.__init__(self, zero)
        basis = {ChainCoset(g, subgroup) for g in group}
        self._basis = list(sorted([b for b in basis if not b == self._zero]))
    def subgroup(self, subgroup):
        return ChainQuotient(subgroup, self._subgroup)
    def __repr__(self):
        return '<%s>' % ', '.join(map(str, self._basis))




if __name__ == '__main__':
    # S = [[0,1,2,3], [3,4],[2,4], [5,6,7]]#, [5,6,7],[6,7,8],[5,7,8],[5,6,8]]
    S = [[0,1,2], [3,4],[4,5],[3,5]]
    # S = [[0,1],[1,2],[0,2],[1,2,3]]
    # S = [[0,1], [2,3]]
    K = SimplicialComplex(lmap(Simplex, S))
    # K = SimplicialComplex(list(map(Simplex, [[0,1,2], [1,2,3]])))
    # L = SimplicialComplex(list(map(Simplex, [[0,1],[0,2],[2,3],[1,3]])))
    C = {d : ChainGroup.simplicial_init(K, d) for d in range(4)}
    # c1 = C[1].subgroup([C[1]._basis[i] for i in (0,1,3,4)])

    # D0 = ChainBoundary(C[0], C[0])
    D1 = ChainBoundary(C[1], C[0])
    D2 = ChainBoundary(C[2], C[1])
    D3 = ChainBoundary(C[3], C[2])

    H0 = ChainQuotient(C[0], D1.im)
    H1 = ChainQuotient(D1.ker, D2.im)

    # print('H0: %d' % len({coset_reduce(c, D1.im) for c in C[0]}))
    # print('H1: %d' % len({coset_reduce(c, D2.im) for c in D1.ker}))
    # print('H2: %d' % len({coset_reduce(c, D3.im) for c in D2.ker}))
