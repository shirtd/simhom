from simhom.algebra import FinitelyGeneratedCoset, FinitelyGeneratedGroup, Quotient
from simhom.util import cref, lmap, lfilt, hstack
from simhom.functions import Homomorphism
from simhom.simplicial import Chain


class ChainGroup(FinitelyGeneratedGroup):
    def __init__(self, K, dim, basis=None):
        self.K, self.dim = K, dim
        self.S = self.K[self.dim]
        if basis is None:
            basis = [Chain({s : 1}) for s in self.S]
        FinitelyGeneratedGroup.__init__(self, basis, Chain(dim=dim))
    def as_reduced(self):
        basis = cref(self.to_mat()).columnspace()
        return self.subgroup(lmap(self.to_element, basis))
    def to_vec(self, c):
        return self.K.to_vec(c)
    def subgroup(self, basis):
        return FinitelyGeneratedGroup(basis, self.zero)
    def to_element(self, v):
        return Chain({self.S[i] : o for i, o in enumerate(v)})
    def subgroup(self, basis=None, f=None):
        basis = self.basis if basis is None else basis
        if f is not None:
            basis = lfilt(f, basis)
        return ChainGroup(self.K, self.dim, basis)
    def __truediv__(self, H):
        assert H.dim == self.dim
        return HomologyGroup(self, H)


class RelativeChainGroup(ChainGroup):
    pass


class Boundary(Homomorphism):
    def __init__(self, X, Y):
        assert (isinstance(X, ChainGroup)
            and isinstance(Y, ChainGroup)
            and Y.dim == X.dim - 1)
        Homomorphism.__init__(self, X, Y)
    def __call__(self, x):
        return x.boundary()
    def __repr__(self):
        return self.im.to_mat().__repr__()


class ChainComplex:
    pass


class RelativeChainComplex(ChainComplex):
    pass


class HomologyClass(FinitelyGeneratedCoset):
    def __init__(self, g, H):
        FinitelyGeneratedCoset.__init__(self, g, H)
    def reduce(self, S):
        E, imap = self.get_imap(S)
        B = cref(hstack(c.as_vec(imap) for c in S))
        return {Chain({E[i] : v for i,v in enumerate(c)}) for c in B.columnspace()}
    def __add__(self, S):
        return HomologyClass(self.g + S.g, self.H)


class HomologyGroup(Quotient):
    def __init__(self, G, H):
        Quotient.__init__(self, G, H)
    def get_coset(self, g):
        return HomologyClass(g, self.H)
