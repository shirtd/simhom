# from simhom.algebra import FGCoset,Quotient, FGZero, FGZeroCoset
from simhom.simplicial import Chain, ZeroChain #, Complex, ComplexElement
from simhom.algebra import FGGroup, FGCoset, FGZeroCoset, Quotient, ZeroElement, ZeroGroup
from simhom.util import cref, lmap, lfilt, hstack, pad
from simhom.functions import Homomorphism, Function
from sympy import zeros


''''''''''''''''''
''' PRIMITIVES '''
''''''''''''''''''


class ChainCoset(FGCoset):
    atom_t = Chain
    def __init__(self, g, H):
        FGCoset.__init__(self, g, H)
    def __invert__(self):
        return ChainCoset(~g, self.H)
    def to_vec(self, g):
        v = zeros(len(self.elems), 1)
        for s,o in g.items():
            v[self.imap[s]] = o
        return v
    def to_element(self, v):
        return Chain({self.elems[i] : o for i, o in enumerate(v)})

class ZeroChainCoset(ChainCoset, FGZeroCoset):
    def __init__(self, G, H):
        ChainCoset.__init__(self, G.zero, H)
        FGZeroCoset.__init__(self, G, H)
    # def __add__(self, g):
    #     return ZeroElement.__add__(self, g)
    # def inverse(self):
    #     return ZeroElement.__add__(self, g)
    # def __eq__(self, g):
    #     return (ZeroElement.__contains__(self, g)
    #         or ChainCoset.__contains__(self, g))
    # def __hash__(self):
    #     return ZeroElement.__hash__(self)


''''''''''''''''''
''' STRUCTURES '''
''''''''''''''''''


class ChainQuotient(Quotient):
    element_t = ChainCoset
    zero_t = ZeroChainCoset
    def __init__(self, G, H):
        self.K, self.dim = G.K, G.dim
        self.S = G.S
        Quotient.__init__(self, G, H)
    def to_vec(self, c):
        return self.K.to_vec(c.g)
    def to_element(self, v):
        c = Chain({self.S[i] : o for i, o in enumerate(v)})
        return ChainCoset(c, self.zero)


class ChainGroup(FGGroup):#, ComplexElement):
    element_t = Chain
    zero_t = ZeroChain
    def __init__(self, K, dim, basis=None, subgroup=None):
        self.K, self.dim = K, dim
        self.S = self.K[self.dim]
        if basis is None:
            basis = [Chain({s : 1}) for s in self.S]
        zero = ZeroChain(subgroup, self.dim)
        FGGroup.__init__(self, basis, zero)
    def as_reduced(self):
        basis = cref(self.to_mat()).columnspace()
        return self.subgroup(lmap(self.to_element, basis))
    def to_vec(self, c):
        return self.K.to_vec(c)
    def subgroup(self, basis):
        return ChainGroup(basis, self.zero)
    def to_element(self, v):
        return Chain({self.S[i] : o for i, o in enumerate(v)})
    def subgroup(self, basis=None, f=None):
        basis = self.basis if basis is None else basis
        if f is not None:
            basis = lfilt(f, basis)
        return ChainGroup(self.K, self.dim, basis, self.zero.subgroup)
    def __truediv__(self, H):
        return ChainQuotient(self, H)
    def zero_group(self):
        return ChainGroup(self.K, self.dim, subgroup=self)
    def as_pair(self, H):
        return ChainGroup(self.K, self.dim, subgroup=H)


# class ZeroChainGroup(ZeroGroup, ChainGroup):
#     def __init__(self, K, dim, subgroup=None):
#         zero = ZeroChain(subgroup, dim)
#         FGGroup.__init__(self, [], zero)


# class ChainComplex: #(Complex):
#     def __init__(self, K, subcomplex=None):
#         self.K, self.subcomplex, self.dim = K, subcomplex, K.dim
#         self._C = {d : ChainGroup(K, d, subgroup=self._subgroup(d)) for d in range(K.dim+2)}
#         self._D = {d : Boundary(self[d], self[d-1]) for d in range(K.dim+2)}
#     def _mkgroup(self, d):
#         return ChainGroup(self.K, d, subgroup=self._subgroup(d))
#     def _subgroup(self, d):
#         if self.subcomplex is not None:
#             return self.subcomplex[d]
#         return None
#     def __getitem__(self, d):
#         if d in self._C:
#             return self._C[d]
#         return self._mkgroup(d)
#     def __call__(self, d):
#         if d in self._D:
#             return self._D[d]
#         return Boundary(self[d], self[d-1])
#     def boundaries(self, d):
#         return self(d+1).im
#     def cycles(self, d):
#         return self(d).ker
#
# class Homology:
#     def __init__(self, C):
#         self.C, self.dim = C, C.dim
#         self._H = {d : C.cycles(d) / C.boundaries(d) for d in range(self.dim+1)}
#     def __getitem__(self, d):
#         if d in self._H:
#             return self._H[d]
#         return self.C.cycles(d) / self.C.boundaries(d)
#     def __repr__(self):
#         pre = {d : 'H%d = span(' % d for d in range(self.dim+1)}
#         return ',\n'.join(['%s%s)' % (pre[d],pad(str(self[d]), len(pre[d])))\
#                             for d in range(self.dim+1)])
