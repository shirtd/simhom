# from simhom.simplicial import Chain
from collections.abc import Hashable, Set
from simhom.util import hstack, stuple, cref
from itertools import combinations
from sympy import Matrix, zeros
from abc import abstractmethod


class GroupElement(Hashable):
    def __hash__(self):
        return hash(self.id)
    @abstractmethod
    def __eq__(self, e):
        pass
    @abstractmethod
    def __add__(self, b):
        pass


class FinitelyGeneratedGroupElement(GroupElement):
    @abstractmethod
    def as_vec(self, imap):
        pass


class Group(Set):
    def __init__(self, basis, zero):
        self.element_t = type(zero)
        if zero in basis:
            basis = [g for g in basis if not g == zero]
        self.elements = set(basis).union({zero})
        self.basis, self.zero = basis, zero
    def is_element(self, e):
        try:
            return isinstance(e, self.element_t) and e in self
        except:
            pass
        return False
    def __iter__(self):
        for v in self.basis:
            yield v
    def __len__(self):
        return len(self.basis)
    def __contains__(self, g):
        return g in self.basis


class FinitelyGeneratedGroup(Group):
    def __init__(self, basis, zero):
        Group.__init__(self, basis, zero)
        assert all(a + b in self for a,b in combinations(self, 2))
    @abstractmethod
    def to_vec(self, e):
        pass
    @abstractmethod
    def to_element(self, v):
        pass
    @abstractmethod
    def subgroup(self, basis):
        pass
    def reduce(self, S):
        M = cref(hstack(self.to_vec(s) for s in S))
        return [self.to_element(c) for c in M.columnspace()]
    def to_mat(self):
        if all(g == self.zero for g in self):
            return zeros(0, len(self))
        return hstack(*[self.to_vec(c) for c in self])
    def __contains__(self, g):
        if not len(self):
            return False
        if isinstance(g, FinitelyGeneratedGroup):
            return all(h in self for h in g)
        M = hstack(self.to_mat(), self.to_vec(g))
        return all(p < len(self) for p in M.rref()[1])
    def __truediv__(self, H):
        return Quotient(self, H)
    def __eq__(self, H):
        return H in self and self in H
    def null(self):
        null = [self.to_element(v) for v in self.to_mat().nullspace()]
        return self.subgroup(null)
    def __repr__(self):
        return 'span(%s)' % ', '.join(map(str, self))


class Coset(Set):
    def __init__(self, g, H):
        self.g, self.H = g, H
        self.elements = {g + h for h in H.elements}
    def __len__(self):
        return len(self.elements)
    def __contains__(self, x):
        return x in self.elements
    def __iter__(self):
        for x in self.elements:
            yield x
    def __eq__(self, Y):
        return all(y in self for y in Y)\
            and all(x in Y for x in self)


class FinitelyGeneratedCoset(Coset, FinitelyGeneratedGroupElement):
    def __init__(self, g, H):
        Coset.__init__(self, g, H.as_reduced())
        self.elements = self.reduce(self.elements) - self.H.elements
        self.id = stuple(c.id for c in self)
    @abstractmethod
    def reduce(self, S):
        pass
    def as_vec(self, imap):
        pass
    def get_imap(self, elements):
        E = stuple({s for c in elements for s in c})
        return E, {s : i for i,s in enumerate(E)}
    def __add__(self, S):
        return FinitelyGeneratedCoset(self.g + S.g, self.H)
    def __eq__(self, S):
        return all(s in S for s in self) and all(s in self for s in S)
    def __contains__(self, c):
        S = list(self.elements) + [c]
        E, imap = self.get_imap(S)
        B, P = hstack(s.as_vec(imap) for s in S).rref()
        return any(p < len(self.elements) for p in P)
    def __hash__(self):
        return hash(self.id)
    def __repr__(self):
        reps = list(sorted(self.elements))
        if not len(reps):
            return '[0]'
        s = '[%s]' % str(reps[0])
        if len(reps) == 1:
            return s
        return '%s = %s' % (s, ('\n%s = ' % (' '*(len(s)))).\
                    join(['[%s]' % str(c) for c in reps[1:]]))


class Quotient(FinitelyGeneratedGroup):
    def __init__(self, G, H):
        self.G, self.H = G, H
        basis, E = set(), set(G.basis)
        zero = self.get_coset(G.zero)
        E -= zero.elements
        while len(E):
            S = self.get_coset(E.pop())
            basis.add(S)
            E -= S.elements
        FinitelyGeneratedGroup.__init__(self, basis, zero)
    @abstractmethod
    def get_coset(self, g):
        pass
    def to_vec(self, e):
        return self.G.to_vec(e.g)
    def to_element(self, v):
        return FinitelyGeneratedCoset(self.G.to_element(v), self.H)
    def subgroup(self, basis):
        return self.G.subgroup(basis) / self.H
    def __contains__(self, S):
        return S.H == self.H and S.g in self.G
    def __repr__(self):
        return ',\n'.join(map(str, self.elements))
