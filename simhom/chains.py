from __future__ import annotations

from simhom.simplicial import *

from simhom.util import stuple, lmap, lfilt, cref, slist
from sympy import Matrix, zeros
from abc import abstractmethod


def csum(*C):
    try:
        dim = next(iter(C)).dim
    except:
        dim = -1
    return sum(*C, ZeroChain(dim=dim))


class AbstractChain(AbstractElement, Atom):
    def __init__(self, id, dim):
        AbstractElement.__init__(self, id)
        self.dim = dim
    def __eq__(self, h) -> bool:
        return (AbstractElement.__eq__(self, h)
            or AbstractElement.__eq__(self, ~h))
    def __hash__(self) -> int:
        return AbstractElement.__hash__(self)
    def __repr__(self) -> str:
        return '%d-Chain' % self.dim
    @abstractmethod
    def __mul__(self, a: int) -> AbstractChain:
        pass


class ZeroChain(ZeroElement, AbstractChain):
    @classmethod
    def is_zero(cls, map={}, dim=-1):
        return not len(map) or all(v == 0 for v in map.values())
    def __init__(self, dim: int=-1, subgroup=None) -> None:
        ZeroElement.__init__(self)
        AbstractChain.__init__(self, (), dim)
    def __len__(self):
        return 0
    def __iter__(self) -> Iterator[Simplex]:
        return iter([])
    def items(self) -> ItemsView[Simplex, int]:
        return {}.items()
    def __repr__(self) -> str:
        return '%s[]' % AbstractChain.__repr__(self)
    def __getitem__(self, s) -> int:
        if isinstance(s, Simplex):
            return 0
        else:
            raise TypeError('ZeroChain contains no simplices')
    def __hash__(self) -> int:
        return AbstractChain.__hash__(self)
    def __mul__(self, a: int) -> AbstractChain:
        return self


class Chain(Generator[Simplex], GroupElement, AbstractChain):
    zero_t = ZeroChain
    @classmethod
    def _mkzero(cls, map={}, dim=-1):
        return ZeroChain(dim)
    @classmethod
    def boundary_init(cls, s: Simplex) -> AbstractChain:
        return Chain({s : (-1)**i for i,s in enumerate(s.face_it())}, s.dim-1)
    def __init__(self, map: Mapping[Simplex, int]={}, dim: int=-1) -> None:
        Generator.__init__(self, map)
        id = tuple((a, self[a] * self[self.get_atom(0)]) for a in self)
        dim = self._atoms[0].dim if len(self) else dim
        assert all(s.dim == dim for s in self)
        AbstractChain.__init__(self, id, dim)
    def __add__(self, c) -> Chain:
        if isinstance(c, ZeroChain):
            return self
        simplices = set(self).union(set(c))
        map = {s : self[s] + c[s] for s in simplices}
        return Chain(map, self.dim)
    def __mul__(self, a: int) -> AbstractChain:
        map = {s : o * a for s,o in self.items()}
        return Chain(map, self.dim)
    def __invert__(self) -> AbstractChain:
        return self * -1
    def normalize(self):
        return self * self[min(self)]
    def __repr__(self) -> str:
        s = ''
        for i, (t, o) in enumerate(self.items()):
            if i > 0 and o > 0:
                s += ' + '
            if o < 0:
                s += ' - '
            if abs(o) > 1:
                s += '%d*' % abs(o)
            s += str(t)
        return '%s[%s]' % (AbstractChain.__repr__(self), s)
    def __lt__(self, other) -> bool:
        return (len(self) < len(other)
            or AbstractElement.__lt__(self, other))
    def __hash__(self) -> int:
        return AbstractChain.__hash__(self)


class ChainGroup(Group[AbstractChain], Basis[Simplex]):
    @classmethod
    def simplicial_init(cls, K, dim):
        basis = [Chain({s : 1}, dim) for s in K[dim]]
        return cls(basis, ZeroChain(dim))
    # @classmethod
    # def reduce(cls, basis):
    #     elems = list({e for b in basis for e in b})
    #     imap = {e : i for i,e in enumerate(elems)}
    #     B = cref(hstack(c.as_vec(imap) for c in basis)).columnspace()
    #     return [Chain({elems[i] : c for i,c in enumerate(v)}) for v in B]
    def __init__(self, basis: AbstractSet[Chain], zero: ZeroChain, dim=-1) -> None:
        self.dim = zero.dim if dim < 0 else dim
        Group.__init__(self, zero)
        Basis.__init__(self, {b for b in basis if b != zero})
    # def to_element(self, v: Matrix) -> Chain:
    #     return Chain({self.get_atom(i) : c for i,c in enumerate(v)})
    def subgroup(self, basis):
        return ChainGroup(basis, self._zero)
    def reduce(self, basis):
        elems = list({e for b in basis for e in b})
        imap = {e : i for i,e in enumerate(elems)}
        B = cref(hstack(c.as_vec(imap) for c in basis)).columnspace()
        return [Chain({elems[i] : c for i,c in enumerate(v)}) for v in B]
    # def reduce(self, S):
    #     C, h = list(sorted(set(S))), 0
    #     for s in self._elems:
    #         if h < len(C):
    #             x = max(range(h, len(C)), key=lambda i: abs(C[i][s]))
    #             if C[x][s] != 0:
    #                 C[h], C[x] = C[x], C[h]
    #                 for i in range(len(C)):
    #                     if i != h:
    #                         C[i] = C[i] - C[h] * (C[i][s] / C[h][s])
    #                 h += 1
    #         else: break
    #     return [c * c[min(c)] for c in sorted(C) if c != 0]


class ChainBoundary(Homomorphism[Chain]):
    element_t, group_t = Chain, ChainGroup
    def __init__(self, X, Y):
        Homomorphism.__init__(self, X, Y)
        self.im = self(self.X)
        self.ker = self._kernel()
    def _element_map(self, x):
        return csum(Chain.boundary_init(s) * o for s, o in x.items())
    def _group_map(self, U):
        return self.Y.subgroup(map(self, U))
    def _kernel(self):
        # if self.im.is_zero():
        #     return self.X
        im = list(map(self, self.X)) + list(self.Y._zero)
        null = hstack(self.Y.to_vec(c) for c in im).nullspace()
        ker = [csum(self.X._basis[i] * v for i,v in enumerate(m[:len(self.X),0])) for m in null]
        return self.X.subgroup(ker)
