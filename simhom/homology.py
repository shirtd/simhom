from __future__ import annotations

from simhom.chains import *


class ChainGroup(Group[AbstractChain], Basis[Simplex]):
    element_t = AbstractChain
    @classmethod
    def simplicial_init(cls, K, dim):
        basis = [Chain({s : 1}, dim) for s in K[dim]]
        return cls(basis, ZeroChain(dim))
    def __init__(self, basis: AbstractSet[Chain], zero: ZeroChain, dim=-1) -> None:
        self.dim = zero.dim if dim < 0 else dim
        Group.__init__(self, zero)
        Basis.__init__(self, {b for b in basis if b != zero})
        self._lmap = {max(c) : c for c in self._basis if len(c)}
        self.depth = 0
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
            else:
                return ChainProjection(self, other)
        if self == 0:
            return ChainInclusion(self, other)
        if other == 0:
            return ChainProjection(self, other)
        raise NotImplementedError('%d-%s >> %d-%s.\n (%s >> %s)' \
                        % (self.dim, str(self), other.dim, str(other),
                            str(type(self)), str(type(other))))
    # def __rrshift__(self, other):
    #     if other == 0:
    #         return ZeroChainInclusion(self)
    #     raise NotImplementedError('%s >> %d-%s.' \
    #         % (str(other), self.dim, str(type(self))))


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
    def __repr__(self):
        return '|bndy>'


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
    def __repr__(self):
        return '|incl>'

class ChainProjection(Homomorphism[AbstractChain]):
    element_t, group_t = AbstractChain, ChainGroup
    def __init__(self, X, Y):
        Homomorphism.__init__(self, X, Y)
        self.im = self(self.X)
        self.ker = self._kernel()
    def _element_map(self, x):
        return x if x in self.Y else self.Y._zero
    def _group_map(self, U):
        return self.Y.subgroup(basis=map(self, U))
    def _kernel(self, zero_group=None):
        ker = [x for x in self.X if self(x) == 0]
        return self.X.subgroup(ker)
    def __repr__(self):
        return '|proj>'

# class ZeroChainInclusion(Homomorphism[AbstractChain]):
#     def __init__(self, Y):
#         self.X, self.Y = 0, Y
#         self.im = Y.subgroup([])
#         self.ker = 0
#     def __repr__(self):
#         return '|incl>'
#
# class ZeroChainProjection(Homomorphism[AbstractChain]):
#     def __init__(self, X):
#         self.X, self.Y = X, 0
#         self.im = 0
#         self.ker = X
#     def _group_map(self, U):
#         return self.Y.subgroup(U)
#     def __repr__(self):
#         return '|proj>'



class ChainQuotient(ChainGroup):
    element_t = AbstractChainCoset
    def __init__(self, group, subgroup, basis=None, zero=None) -> None:
        self._group, self._subgroup, self.dim = group, subgroup, subgroup.dim
        zero = ZeroChainCoset(group._zero, subgroup) if zero is None else zero
        basis = {ChainCoset(g, subgroup) for g in group} if basis is None else basis
        ChainGroup.__init__(self, [c for c in basis if not c == zero], zero)
        self.depth = self._subgroup.depth + 1
    def subgroup(self, group=None, basis=None):
        group = self._group if group is None else group
        # return ChainQuotient(basis, self._zero, self._subgroup)
        return ChainQuotient(group, self._subgroup, basis, self._zero)
    def reduce(self, basis):
        return list(sorted(basis))
    def __contains__(self, c):
        return (any(c == l for l in self)
            or c == self._zero or c == 0)
    def __rshift__(self, other):
        if other == 0:
            try:
                zero_group = ChainGroup([], self._group._zero)
                return InducedQuotientMap(self, zero_group / self._subgroup)
            except AttributeError as e:
                print(self, type(self), self.dim, self._group)
                raise e
        elif self.depth == 1 and other.depth == 0:
            return ConnectingHomomorphism(self, other)
        return InducedQuotientMap(self, other)

class ConnectingHomomorphism(Homomorphism[Chain]):
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
        basis = list(map(self, self.X._group))
        if zero_group is not None:
            basis += list(zero_group._basis)
        null = to_mat(basis).nullspace()
        ker = [csum(self.X._group._basis[i] * v for i,v in enumerate(m[:len(self.X._group),0])) for m in null]
        return self.X.subgroup(self.X._group.subgroup(ker))
    def __repr__(self):
        return '|conn>'


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
    def __repr__(self):
        return '|proj>'

class ChainQuotientInclusion(Homomorphism[AbstractChain]):
    element_t, group_t = AbstractChain, ChainGroup
    def __init__(self, X, Y):
        Homomorphism.__init__(self, X, Y)
        self.im = self(self.X)
        self.ker = self._kernel()
    def _element_map(self, x):
        return x._rep
    def _group_map(self, U):
        return self.Y.subgroup(map(self, U))
    def _kernel(self, zero_group=None):
        ker = []
        if zero_group is not None:
            ker = [x for x in self.X if x in zero_group]
        return self.X.subgroup(self.X._group(ker))
    def __repr__(self):
        return '|proj>'


class InducedQuotientMap(Homomorphism[AbstractChainCoset]):
    element_t, group_t = AbstractChainCoset, ChainQuotient
    def __init__(self, X, Y):
        Homomorphism.__init__(self, X, Y)
        self._f = X._group >> Y._group
        self.im = self(self.X)
        self.ker = self._kernel()
    def _element_map(self, x):
        return ChainCoset(self._f(x), self.Y._subgroup)
    def _group_map(self, U):
        if isinstance(U, ChainQuotient):
            return self.Y.subgroup(self._f(U._group))
        return self.Y.subgroup(self._f(U))
    def _kernel(self, zero_group=None):
        ker = self._f._kernel(self.Y._zero._subgroup)
        return self.X.subgroup(ker)
    def __repr__(self):
        return '%s*' % self._f.__repr__()
