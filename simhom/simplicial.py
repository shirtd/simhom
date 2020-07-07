
from __future__ import annotations

from simhom.algebra import *

# from simhom.algebra import BasisElement, ZeroElement
from simhom.util import stuple, lmap, lfilt, slist
# from collections.abc import Mapping, MutableMapping
from typing import List
# from sympy import Matrix, zeros
from abc import abstractmethod


class Simplex(tuple, Atom):
    def __new__(cls, *v: Union[int, Collection[int]]) -> Simplex:
        try:
            return tuple.__new__(cls, stuple(*iter(v)))
        except:
            return tuple.__new__(cls, stuple(v))
    def __init__(self, *v: Union[int, Collection[int]]) -> None:
        tuple.__init__(self)
        self.dim = len(self) - 1
    def face_it(self) -> Iterator[Simplex]:
        for i in range(self.dim+1):
            yield Simplex(self[:i]+self[i+1:])
    def __contains__(self, s) -> bool:
        if isinstance(s, Simplex):
            return all(v in self for v in s)
        elif isinstance(s, int):
            return tuple.__contains__(self, s)
        return False


class SimplicialSkeleton:
    def __init__(self, dim: int=-1) -> None:
        self.S: List[Simplex] = []
        self.imap: Dict[Simplex, int] = {}
        self.dim = dim
    def __len__(self) -> int:
        return len(self.S)
    def __getitem__(self, s) -> Union[Simplex, int]:
        if isinstance(s, Simplex):
            return self.imap[s]
        return self.S[s]
    def __iter__(self) -> Iterator[Simplex]:
        for s in self.S:
            yield s
    def add(self, s: Simplex) -> None:
        if self.dim == -1:
            self.dim = s.dim
        elif self.dim != s.dim:
            raise Exception('dimension of %s does not match' % str(s))
        if not s in self:
            self.imap[s] = len(self)
            self.S.append(s)
    def __contains__(self, c: Simplex) -> bool:
        return c in self.imap
    def __repr__(self) -> str:
        return str(self.S)

''''''''''''''''''
''' STRUCTURES '''
''''''''''''''''''


class SimplicialComplex:
    def __init__(self, S: List[List[int]]) -> None:
        self.dim = max(len(s) for s in S) - 1
        self.S = {d : SimplicialSkeleton(d) for d in range(self.dim+1)}
        for s in S:
            self.add(Simplex(s))
    def add(self, s: Simplex) -> None:
        self[s.dim].add(s)
        if s.dim > 0:
            for t in s.face_it():
                self.add(t)
    def __getitem__(self, d: int) -> SimplicialSkeleton:
        if d > self.dim or d < 0:
            return SimplicialSkeleton(d)
        return self.S[d]
    def items(self) -> ItemsView[int, SimplicialSkeleton]:
        return self.S.items()
    def __iter__(self) -> Iterator[SimplicialSkeleton]:
        for k,v in self.S.items():
            yield v
    def __contains__(self, s: Simplex) -> bool:
        return s in self[s.dim]
    def __repr__(self) -> str:
        return '{%s}' % ',\n '.join(['%d: %s' % (d, str(s)) for d,s in self.items()])
