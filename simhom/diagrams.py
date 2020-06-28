from __future__ import annotations

from typing import Type, TypeVar, Generic, Callable, Collection, Tuple,\
                    Container, NewType, Hashable, AbstractSet, Dict,\
                     Iterable, Mapping, Union, ItemsView, Iterator
from simhom.util import hstack, stuple, cref, slist
from itertools import combinations
from sympy import Matrix, zeros
from abc import abstractmethod


O = TypeVar('O', bound='Object')
M = TypeVar('M', bound='Morphism')
E = TypeVar('E', bound='GroupElement')
A = TypeVar('A', bound='Atom')

''''''''''''''''''
''' CATEGORIES '''
''''''''''''''''''
'''' OBJECTS '''''
class Object:
    def __init__(self) -> None:
        pass

''''''''''''''''''
'''' MORPHISMS '''
class Morphism(Generic[O]): # Object O
    @classmethod
    def create(cls: Type[M], X: O, Y: O) -> Morphism:
        return cls(X, Y)
    def __init__(self, X : O, Y : O) -> None:
        self.X, self.Y = X, Y
    def __call__(self, x: O) -> O:
        pass


''''''''''''''''''
''''' GROUPS '''''
''''''''''''''''''
''''' ELEMENTS '''
class GroupElement:
    @classmethod
    def create(cls: Type[E], id: tuple) -> GroupElement:
        return cls(id)
    def __init__(self, id: tuple) -> None:
        self.id = id
        self.is_zero = False
    def __eq__(self, other) -> bool:
        return self.id == other.id
    def __hash__(self):
        return hash(self.id)
    def __sub__(self, h):
        return self + ~h
    @abstractmethod
    def __invert__(self):
        pass
    @abstractmethod
    def __add__(self, other):
        pass


class ZeroElement(GroupElement):
    # @classmethod
    # def create(cls: Type[E]) -> ZeroElement:
    #     return cls()
    def __init__(self) -> None:
        GroupElement.__init__(self, tuple())
        self.is_zero = True
    def __add__(self, other : GroupElement) -> GroupElement:
        return other
    def __invert__(self) -> GroupElement:
        return self
    def __eq__(self, other) -> bool:
        return other.is_zero


class Atom:
    pass
    # @classmethod
    # def create(cls: Type[A]) -> Atom:
    #     return cls()


class Generator(Mapping[A, int]): # Atom A
    # @classmethod
    # def create(cls, map: Mapping[A, int]) -> Generator:
    #     return cls(map)
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
class Group(Object, AbstractSet[E]):
    # @classmethod
    # def create(cls: Type[O], zero: E) -> Generator:
    #     return cls(zero)
    def __init__(self, zero : E) -> None:
        Object.__init__(self)
        self.zero = zero

# Finite Group
class FiniteGroup(Group[E]):
    # @classmethod
    # def create(cls: Type[O], zero : E, elements : AbstractSet[E]) -> FiniteGroup:
    #     return cls(zero, elements)
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
    # @classmethod
    # def create(cls: Type[O], basis : AbstractSet[Generator[A]]) -> Basis:
    #     return cls(basis)
    def __init__(self, basis : AbstractSet[Generator[A]]) -> None:
        self._elems = list({e for b in basis for e in b})
        self._imap = {e : i for i,e in enumerate(self._elems)}
        self._basis = basis
    def __len__(self) -> int:
        return len(self._basis)
    def __contains__(self, g) -> bool:
        if not (len(self) or all(e in self._imap for e in g)):
            return False
        M = hstack(self.to_mat(), self.to_vec(g))
        return all(p < len(self._basis) for p in M.rref()[1])
    def to_mat(self) -> Matrix:
        return hstack(map(self.to_vec, self._basis))
    def reduce(self, S: AbstractSet[Generator[A]]) -> AbstractSet[Generator[A]]:
        B = cref(hstack(self.to_vec(c) for c in S))
        return {self.to_element(v) for v in B.columnspace()}
    def __eq__(self, S) -> bool:
        return (all(s in self for s in S)
            and all(s in S for s in self))
    def to_vec(self, g: Generator[A]) -> Matrix:
        return g.as_vec(self._imap)
    def __repr__(self) -> str:
        return '<%s>' % ', '.join(map(str, self))
    @abstractmethod
    def to_element(self, v: Matrix) -> Generator[A]:
        pass

''''''''''''''''''
'''' MORPHISMS '''
class Homomorphism(Morphism[Group[E]]):
    # @classmethod
    # def create(cls: Type[M], X: Group[E], Y: Group[E]) -> Homomorphism:
    #     return cls(X, Y)
    def __init__(self, X: Group[E], Y: Group[E]):
        Morphism.__init__(self, X, Y)


# '''''''''''''''''''''''
# '''' FINITE GROUPS ''''
# '''''''''''''''''''''''
# ''''''' OBJECTS '''''''
#
# class SmallGroupElement(GroupElement):
#     def __init__(self, i : int) -> None:
#         self._i = i % 2
#         GroupElement.__init__(self, (i,))
#     def __add__(self, other):
#         return SmallGroupElement(self._i + other._i)
#     def __repr__(self):
#         return str(self._i)
#
# class SmallGroupZeroElement(ZeroElement, SmallGroupElement):
#     def __init__(self):
#         SmallGroupElement.__init__(self, 0)
#         ZeroElement.__init__(self)
#
#
# class SmallGroup(FiniteGroup[SmallGroupElement]):
#     def __init__(self):
#         zero = SmallGroupZeroElement()
#         elements = {SmallGroupElement(1)}
#         FiniteGroup.__init__(self, zero, elements)
