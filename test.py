from __future__ import annotations

from simhom.diagrams import *


def csum(*C):
    try:
        dim = next(iter(C)).dim
    except:
        dim = -1
    return sum(*C, ZeroChain(dim=dim))


class Simplex(tuple, Atom):
    def __new__(cls, *v: Union[int, Collection[int]]) -> Simplex:
        try:
            return tuple.__new__(cls, stuple(*iter(v)))
        except:
            return tuple.__new__(cls, stuple(v))
    def __init__(self, *v: Union[int, Collection[int]]) -> None:
        tuple.__init__(self)
        self.dim = len(self) - 1
    def boundary(self) -> Chain:
        c = {}
        if self.dim > 0:
            c = {self.get_face(i) : (-1)**i for i in range(self.dim+1)}
        return Chain(c, self.dim-1)
    def get_face(self, i: int) -> Simplex:
        return Simplex(self[:i]+self[i+1:])
    def __contains__(self, s) -> bool:
        if isinstance(s, Simplex) or isinstance(s, Chain):
            return all(v in self for v in s)
        return s in self


class AbstractChain(GroupElement):
    def __init__(self, id, dim):
        GroupElement.__init__(self, id)
        self.dim = dim
    def __eq__(self, h) -> bool:
        return (GroupElement.__eq__(self, h)
            or GroupElement.__eq__(self, ~h))
    def __hash__(self) -> int:
        return GroupElement.__hash__(self)
    def __repr__(self) -> str:
        return '%d-Chain' % self.dim
    @abstractmethod
    def boundary(self) -> AbstractChain:
        pass


class Chain(Generator[Simplex], AbstractChain):
    def __new__(cls, map: Mapping[Simplex, int]={}, dim: int=-1):
        if not len(map) or all(v == 0 for v in map.values()):
            return ZeroChain(dim)
        return super().__new__(cls)
    def __init__(self, map: Mapping[Simplex, int]={}, dim: int=-1) -> None:
        Generator.__init__(self, map)
        id = tuple((a, self[a] * self[self.get_atom(0)]) for a in self)
        dim = self._atoms[0].dim if len(self) else dim
        assert all(s.dim == dim for s in self)
        AbstractChain.__init__(self, id, dim)
    def __add__(self, c) -> Chain:
        simplices = set(self).union(set(c))
        map = {s : self[s] + c[s] for s in simplices}
        return Chain(map, self.dim)
    def __mul__(self, a: int) -> Chain:
        map = {s : o * a for s,o in self.items()}
        return Chain(map, self.dim)
    def __invert__(self) -> Chain:
        return self * -1
    def boundary(self) -> Chain:
        return csum(s.boundary() * o for s, o in self.items())
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


class ZeroChain(ZeroElement, AbstractChain):
    def __init__(self, dim: int=-1, subgroup=None) -> None:
        print('ZeroChain.__init__(self, %s, %s)' % (str(dim), str(subgroup)))
        ZeroElement.__init__(self)
        AbstractChain.__init__(self, (), dim)
    def __iter__(self) -> Iterator[Simplex]:
        return iter([])
    def items(self) -> ItemsView[Simplex, int]:
        return {}.items()
    def boundary(self) -> AbstractChain:
        return ZeroChain(dim=self.dim-1)
    def __repr__(self) -> str:
        return '%s[]' % AbstractChain.__repr__(self)


if __name__ == '__main__':
    s = Simplex(0, 1, 2)
