# from simhom.simplicial import Chain
from collections.abc import Hashable, Set, Mapping
from simhom.util import hstack, stuple, cref, slist
from itertools import combinations
from sympy import Matrix, zeros
from abc import abstractmethod


''''''''''''''''''
''' PRIMITIVES '''
''''''''''''''''''


class GroupElement(Hashable):
    def __init__(self, id):
        self.id = id
    def is_zero(self):
        return isinstance(self, ZeroElement)
    def __hash__(self):
        return hash(self.id)
    def __lt__(self, c):
        return self.id < c.id
    def __eq__(self, h):
        return self.id == h.id
    @abstractmethod
    def inverse(self):
        pass
    @abstractmethod
    def __add__(self, g):
        pass


class ZeroElement(GroupElement):
    def __init__(self, subgroup=None):
        self.subgroup = subgroup
        self.basis = subgroup.basis if subgroup is not None else []
        GroupElement.__init__(self, tuple())
    def __add__(self, g):
        return g
    def inverse(self):
        return self
    def __eq__(self, g):
        return (g.is_zero()
            or (self.subgroup is not None
            and g in self.subgroup))
    def __hash__(self):
        return GroupElement.__hash__(self)


'''
an element that can be vectorized.
    maps atoms to coefficients '''
class BasisElement(GroupElement, Mapping):
    def __init__(self, map):
        if not isinstance(map, dict):
            self.map = {k : 1 for k in map}
        else:
            self.map = {k : v for k,v in map.items() if v != 0}
        self.atoms = stuple(self.map)
        id = tuple((a, self[a]) for a in self)
        GroupElement.__init__(self, id)
    def items(self):
        for k,v in self.map.items():
            yield (k,v)
    def __len__(self):
        return len(self.atoms)
    def __iter__(self):
        for s in self.atoms:
            yield s
    def __getitem__(self, a):
        return self.map[a] if a in self else 0
    def __contains__(self, a):
        return a in self.map
    def __hash__(self):
        return GroupElement.__hash__(self)
    def __sub__(self, h):
        return self + h.inverse()
    def as_vec(self, imap):
        v = zeros(len(imap),1)
        for s,o in self.items():
            v[imap[s]] = o
        return v
    @property
    @abstractmethod
    def atom_t(self):
        pass


''''''''''''''''''
''' STRUCTURES '''
''''''''''''''''''


class Group:
    def __init__(self, zero):
        self.zero = zero
    @abstractmethod
    def __contains__(self, g):
        pass
    @property
    @abstractmethod
    def zero_t(self):
        pass
    @property
    @abstractmethod
    def element_t(self):
        pass


class Basis:
    def __init__(self, basis):
        self.elems = list({e for b in basis for e in b})
        self.imap = {e : i for i,e in enumerate(self.elems)}
        self.basis = basis
    def __contains__(self, g):
        if (not len(self.basis)
            or any(not e in self.imap for e in g)):
            return False
        M = hstack(self.to_mat(), self.to_vec(g))
        return all(p < len(self.basis) for p in M.rref()[1])
    def to_mat(self):
        return hstack(map(self.to_vec, self.basis))
    def reduce(self, S):
        B = cref(hstack(self.to_vec(c) for c in S))
        return {self.to_element(v) for v in B.columnspace()}
    @abstractmethod
    def to_vec(self, g):
        pass
    @abstractmethod
    def to_element(self, v):
        pass


class FGGroup(Group, Basis):
    def __init__(self, basis, zero):
        Group.__init__(self, zero)
        Basis.__init__(self, basis)
    def __len__(self):
        return len(self.basis)
    def __iter__(self):
        for g in self.basis:
            yield g
    def __repr__(self):
        return 'span(%s)' % ', '.join(map(str, self))


class Quotient(FGGroup):
    def __init__(self, G, H):
        basis, zero = set(), self.zero_t(G, H)
        E = set(G.basis) - set(zero.basis)
        while len(E):
            S = self.element_t(E.pop(), zero)
            if zero != S:
                basis.add(S)
            E -= set(S.basis)
        FGGroup.__init__(self, list(basis), zero)
    def __repr__(self):
        classes = [self.zero] + self.basis
        return ',\n'.join(map(str, classes))


''''''''''''''''''
''''' HYBRID '''''
''''''''''''''''''


class FGCoset(Basis, BasisElement):
    def __init__(self, g, H):
        self.g, self.H = g, H
        # Hbasis = H.subgroup.reduce(H.basis)
        elements = {g + h for h in H.basis}.union({g})
        Basis.__init__(self, elements)
        self.basis = self.reduce(self.basis) - self.reduce(H.basis)
        BasisElement.__init__(self, self.basis)
    def __contains__(self, g):
        return Basis.__contains__(self, g)
    def __eq__(self, S):
        return (all(s in self for s in S)
            and all(s in S for s in self))
    def __hash__(self):
        return BasisElement.__hash__(self)
    def __add__(self, S):
        return self.__class__(self.g + S.g, self.H)
    def _mkstr(self, reps):
        if not len(reps):
            return '[0]'
        s = '[%s]' % str(reps[0])
        if len(reps) == 1:
            return s
        return '%s = %s' % (s, ('\n%s = ' % (' '*(len(s)))).\
                    join(['[%s]' % str(c) for c in reps[1:]]))
    def __repr__(self):
        return self._mkstr(slist(self.basis))


class FGZeroCoset(FGCoset, ZeroElement):
    def __init__(self, G, H):
        FGCoset.__init__(self, G.zero, H)
        ZeroElement.__init__(self, H)
    def __repr__(self):
        return self._mkstr(slist([self.g]+self.basis))
