from __future__ import annotations

from simhom.chains import *


class ChainGroup(Group[AbstractChain], Basis[Simplex]):
    @classmethod
    def simplicial_init(cls, K, dim):
        basis = [Chain({s : 1}, dim) for s in K[dim]]
        return cls(basis, ZeroChain(dim))
    def __init__(self, basis: AbstractSet[Chain], zero: ZeroChain, dim=-1) -> None:
        self.dim = zero.dim if dim < 0 else dim
        Group.__init__(self, zero)
        Basis.__init__(self, {b for b in basis if b != zero})
        self._lmap = {max(c) : c for c in self._basis if len(c)}
    def subgroup(self, basis):
        return ChainGroup(basis, self._zero)
    def is_subgroup(self, other):
        return all(c in other for c in self)
    def canonical(self, basis):
        C = list(sorted(basis))
        for h in range(len(C)):
            s = max(C[h])
            C[h] = C[h] * C[h][s]
            for i in range(len(C)):
                if i != h and s in C[i]:
                    C[i] = C[i] - C[h] * (C[i][s] / C[h][s])
        return C
    def reduce(self, S):
        C, h = list(sorted(set(S))), 0
        for s in self._elems:
            if h < len(C):
                x = max(range(h, len(C)), key=lambda i: abs(C[i][s]))
                if C[x][s] != 0:
                    C[h], C[x] = C[x], C[h]
                    for i in range(len(C)):
                        if i != h:
                            C[i] = C[i] - C[h] * (C[i][s] / C[h][s])
                    h += 1
            else: break
        return self.canonical([c * c[min(c)] for c in sorted(C) if c != 0])
    def __truediv__(self, other):
        return ChainQuotient(self, other)
    def __rshift__(self, other):
        if other.dim == self.dim - 1:
            if not isinstance(other, ChainQuotient):
                return ChainBoundary(self, other)
        elif other.dim == self.dim:
            if isinstance(other, ChainQuotient):
                return ChainQuotientProjection(self, other)
            elif self.is_subgroup(other):
                return ChainInclusion(self, other)
        raise NotImplementedError('%d-%s >> %d-%s.' \
                        % (self.dim, str(type(self)),
                            other.dim, str(type(other))))


class ChainBoundary(Homomorphism[Chain]):
    element_t, group_t = AbstractChain, ChainGroup
    def __init__(self, X, Y):
        Homomorphism.__init__(self, X, Y)
        self.im = self(self.X)
        self.ker = self._kernel()
    def _element_map(self, x):
        return csum(Chain.boundary_init(s) * o for s, o in x.items())
    def _group_map(self, U):
        return self.Y.subgroup(map(self, U))
    def _kernel(self, zero_group=None):
        basis = list(map(self, self.X))
        if zero_group is not None:
            basis += list(zero_group._basis)
        null = to_mat(basis).nullspace()
        ker = [csum(self.X._basis[i] * v for i,v in enumerate(m[:len(self.X),0])) for m in null]
        return self.X.subgroup(ker)


class ChainInclusion(Homomorphism[AbstractChain]):
    element_t, group_t = AbstractChain, ChainGroup
    def __init__(self, X, Y):
        Homomorphism.__init__(self, X, Y)
        self.im = self(self.X)
        self.ker = self._kernel()
    def _element_map(self, x):
        return x
    def _group_map(self, U):
        return self.Y.subgroup(U)
    def _kernel(self, zero_group=None):
        ker = []
        if zero_group is not None:
            ker = [x for x in self.X if x in zero_group]
        return self.X.subgroup(ker)


class ChainQuotient(ChainGroup):
    def __init__(self, group, subgroup) -> None:
        self._group, self._subgroup, self.dim = group, subgroup, group.dim
        zero = ZeroChainCoset(group._zero, subgroup)
        basis = {ChainCoset(g, subgroup) for g in group}
        ChainGroup.__init__(self, [c for c in basis if not c == zero], zero)
    def subgroup(self, subgroup):
        return ChainQuotient(subgroup, self._subgroup)
    def reduce(self, basis):
        return list(sorted(basis))
    def __contains__(self, c):
        return (any(c == l for l in self)
            or c == self._zero)
    def __rshift__(self, other):
        return InducedMap(self, other)


class ChainQuotientProjection(Homomorphism[AbstractChain]):
    element_t, group_t = AbstractChain, ChainGroup
    def __init__(self, X, Y):
        Homomorphism.__init__(self, X, Y)
        self.im = self(self.X)
        self.ker = self._kernel()
    def _element_map(self, x):
        return ChainCoset(x, self.Y._subgroup)
    def _group_map(self, U):
        return self.Y.subgroup(U)
    def _kernel(self, zero_group=None):
        ker = [x for x in self.X if self(x) == self.Y._zero]
        return self.X.subgroup(ker)


class InducedMap(Homomorphism[AbstractChainCoset]):
    element_t, group_t = AbstractChainCoset, ChainQuotient
    def __init__(self, X, Y):
        Homomorphism.__init__(self, X, Y)
        self._f = X._group >> Y._group
        self.im = self(self.X)
        self.ker = self._kernel()
    def _element_map(self, x):
        return ChainCoset(self._f(x), self.Y._subgroup)
    def _group_map(self, U):
        return self.Y.subgroup(self._f(self.X._group))
    def _kernel(self):
        ker = self._f._kernel(self.Y._zero._subgroup)
        return self.X.subgroup(ker)
