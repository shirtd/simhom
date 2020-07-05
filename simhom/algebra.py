from __future__ import annotations

from typing import Type, TypeVar, Generic, Callable, Collection, Tuple,\
                    Container, NewType, Hashable, AbstractSet, Dict,\
                     Iterable, Mapping, Union, ItemsView, Iterator
from simhom.util import hstack, stuple, cref, slist
from itertools import combinations
from sympy import Matrix, zeros
from abc import abstractmethod
from simhom.categories import *


# O = TypeVar('O', bound='Object')
# M = TypeVar('M', bound='Morphism')
E = TypeVar('E', bound='AbstractElement')
A = TypeVar('A', bound='Atom')
N = TypeVar('N', bound='Generator')


class AbstractElement:
    def __init__(self, id: tuple) -> None:
        self.id = id
    def __eq__(self, other) -> bool:
        return self.id == other.id
    def __hash__(self):
        return hash(self.id)
    def __sub__(self, h):
        return self + ~h
    def __lt__(self, other):
        return self.id < other.id
    @abstractmethod
    def __invert__(self):
        pass
    @abstractmethod
    def __add__(self, other):
        pass


class ZeroElement(AbstractElement):
    def __init__(self) -> None:
        AbstractElement.__init__(self, tuple())
    def __add__(self, other) -> AbstractElement:
        return other
    def __invert__(self) -> AbstractElement:
        return self
    def __eq__(self, other) -> bool:
        return isinstance(other, ZeroElement) or other == 0
    @classmethod
    @abstractmethod
    def is_zero(cls, *args, **kwargs):
        pass

class GroupElement(AbstractElement):
    def __new__(cls, *args, **kwargs):
        if cls.zero_t.is_zero(*args, **kwargs):
            return cls._mkzero(*args, **kwargs)
        return super().__new__(cls)
    def __init__(self, id):
        AbstractElement.__init__(self, id)
    @property
    @abstractmethod
    def zero_t(self):
        pass
    @abstractmethod
    def _mkzero(self, *args, **kwargs):
        pass


class Atom:
    pass


class Generator(Mapping[A, int]): # Atom A
    def __init__(self, map: Mapping[A, int]):
        self._map = {k : v for k,v in map.items() if v != 0}
        self._atoms = stuple(self._map)
    def items(self) -> ItemsView[A, int]:
        return self._map.items()
    def __len__(self) -> int:
        return len(self._atoms)
    def __iter__(self) -> Iterator[A]:
        for s in self._atoms:
            yield s
    def __getitem__(self, a: A) -> int:
        return self._map[a] if a in self else 0
    def get_atom(self, i: int) -> A:
        return self._atoms[i]
    def __contains__(self, a) -> bool:
        return a in self._map
    def as_vec(self, imap : Mapping[A, int]) -> Matrix:
        v = zeros(len(imap), 1)
        for s, o in self.items():
            v[imap[s]] = o
        return v


''''''''''''''''
'''' OBJECTS '''
class Group(Object, Generic[E]):
    def __init__(self, zero : E) -> None:
        Object.__init__(self)
        self._zero = zero
    @abstractmethod
    def subgroup(self, *args, **kwargs):
        pass

# Finite Group
class FiniteGroup(Group[E]):
    def __init__(self, zero : E, elements : AbstractSet[E]) -> None:
        Group.__init__(self, zero)
        self._elements = {zero}.union(elements)
    def __contains__(self, e) -> bool:
        return e in self._elements
    def __iter__(self) -> Iterator[E]:
        for e in self._elements:
            yield e
    def __len__(self) -> int:
        return len(self._elements)

class Basis(AbstractSet[Generator[A]]):
    def __init__(self, basis : AbstractSet[Generator[A]]) -> None:
        self._elems = list({e for b in basis for e in b})
        self._imap = {e : i for i,e in enumerate(self._elems)}
        self._basis = self.reduce(basis)
    def __len__(self) -> int:
        return len(self._basis)
    def __eq__(self, S) -> bool:
        return (all(s in self for s in S)
            and all(s in S for s in self))
    def __iter__(self) -> Iterator:
        for b in self._basis:
            yield b
    def get_atom(self, i):
        return self._elems[i]
    def has_atom(self, a):
        return a in self._elems
    def to_vec(self, g: Generator[A]) -> Matrix:
        return g.as_vec(self._imap)
    def __repr__(self) -> str:
        return '<%s>' % ', '.join(map(str, self))
    def __contains__(self, other):
        if not all(self.has_atom(a) for a in other):
            return False
        R = self.reduce(self._basis + [other])
        return len(R) <= len(self)
    @abstractmethod
    def reduce(self, S):
        pass
    # @abstractmethod
    # def to_element(self, v: Matrix) -> Generator[A]:
    #     pass
    # def __contains__(self, g) -> bool:
    #     if not (len(self) and all(e in self._imap for e in g)):
    #         return False
    #     M = hstack(self.to_mat(), self.to_vec(g))
    #     return all(p < len(self._basis) for p in M.rref()[1])
    # def to_mat(self) -> Matrix:
    #     return hstack(map(self.to_vec, self._basis))
    # def reduce(self, S: AbstractSet[Generator[A]]) -> AbstractSet[Generator[A]]:
    #     B = cref(hstack(self.to_vec(c) for c in S))
    #     return {self.to_element(v) for v in B.columnspace()}

''''''''''''''''''
'''' MORPHISMS '''
class Homomorphism(Morphism[Group[E]]):
    def __init__(self, X: Group[E], Y: Group[E]) -> None:
        Morphism.__init__(self, X, Y)
    def __call__(self, x: Union[AbstractElement, Group[E]])\
            -> Union[AbstractElement, Group[E]]:
        if isinstance(x, AbstractElement):
            return self._element_map(x)
        elif isinstance(x, Group):
            return self._group_map(x)
        else:
            raise TypeError(str(x))
    @abstractmethod
    def _element_map(self, x: AbstractElement) -> AbstractElement:
        pass
    @abstractmethod
    def _group_map(self, U: Group[E]) -> Group[E]:
        pass
