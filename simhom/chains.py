from __future__ import annotations

from simhom.simplicial import *

from simhom.util import stuple, lmap, lfilt, slist, to_mat
# from sympy import Matrix, zeros
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
    # def __repr__(self) -> str:
    #     return '%d-Chain' % self.dim
    @abstractmethod
    def __mul__(self, a: int) -> AbstractChain:
        pass


class ZeroChain(ZeroElement, AbstractChain):
    @classmethod
    def is_zero(cls, map={}, dim=-1):
        return not len(map) or all(v == 0 for v in map.values())
    def __init__(self, dim: int=-1) -> None:
        ZeroElement.__init__(self)
        AbstractChain.__init__(self, (), dim)
        self._map = {}
    def __len__(self):
        return 0
    def __iter__(self) -> Iterator[Simplex]:
        return iter([])
    def items(self) -> ItemsView[Simplex, int]:
        return {}.items()
    def __repr__(self) -> str:
        # return '%s[]' % AbstractChain.__repr__(self)
        return 'O'
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
        # return '%s[%s]' % (AbstractChain.__repr__(self), s)
        return '%s' % s
    def __lt__(self, other) -> bool:
        return (len(self) < len(other)
            or AbstractElement.__lt__(self, other))
    def __hash__(self) -> int:
        return AbstractChain.__hash__(self)


class AbstractChainCoset(AbstractChain):
    def __init__(self, chain, subgroup):
        for s in chain:
            if s in subgroup._lmap:
                c = subgroup._lmap[s]
                chain = chain - c * (chain[s] / c[s])
        self._rep = chain
        self._subgroup = subgroup
        AbstractChain.__init__(self, self._rep.id, self._rep.dim)
    def __eq__(self, other):
        try:
            return other - self._rep in self._subgroup
        except TypeError as e:
            print(other, self._rep)
            print(self._subgroup)
            raise e
    def __hash__(self):
        return AbstractChain.__hash__(self)


class ZeroChainCoset(AbstractChainCoset, ZeroChain):
    @classmethod
    def is_zero(cls, chain, subgroup):
        return chain in subgroup
    def __init__(self, chain, subgroup):
        AbstractChainCoset.__init__(self, chain, subgroup)
        ZeroChain.__init__(self, chain.dim)
    def __eq__(self, other):
        return (ZeroChain.__eq__(self, other)
            or AbstractChainCoset.__eq__(self, other))
    def __hash__(self):
        return AbstractChainCoset.__hash__(self)
    def __mul__(self, a: int) -> AbstractChain:
        return self
    def __repr__(self) -> str:
        return '[O]' % AbstractChain.__repr__(self)
    def __iter__(self):
        return iter(self._rep)


class ChainCoset(AbstractChainCoset, Chain):
    zero_t = ZeroChainCoset
    @classmethod
    def _mkzero(cls, chain, subgroup):
        return ZeroChainCoset(chain, subgroup)
    def __init__(self, chain, subgroup):
        AbstractChainCoset.__init__(self, chain, subgroup)
        Chain.__init__(self, self._rep._map, self._rep.dim)
    def __add__(self, other):
        return ChainCoset(Chain.__add__(self, other), self._subgroup)
    def __mul__(self, a):
        return ChainCoset(self._rep * a, self._subgroup)
    def __repr__(self):
        return '[%s]' % str(self._rep)
