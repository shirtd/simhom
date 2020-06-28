# from simhom.algebra import FGCoset,Quotient, FGZero, FGZeroCoset
from simhom.simplicial import Chain, ZeroChain #, Complex, ComplexElement
from simhom.algebra import FGGroup, FGCoset, FGZeroCoset, Quotient, ZeroElement
from simhom.homology import ChainGroup, ChainQuotient
from simhom.util import cref, lmap, lfilt, hstack, pad
from simhom.functions import Identity, Homomorphism, Epimorphism, Monomorphism, Function
from collections.abc import Sequence
from abc import abstractmethod
from sympy import zeros


# class ExactSequenceElement:
#     def __init__(self, next):
#         self.next = ExactSequenceElement(it)
#
# class ChainGroupElement(ChainGroup, ExactSequenceElement):
#     def __init__(self, next, *args, **kwargs):
#         ChainGroup.__init__(self, *args, **kwargs)
#         ExactSequenceElement.__init__(self, next)
#
# class ExactSequence:
#     def __init__(self, sequence):
#         last = None
#         self._sequence = [last]
#         for args in reversed(sequence):
#             last = self._mkelem(args)
#             self._sequence.insert(0, last)
#         self._sequence.insert(0, None)
#     @abstractmethod
#     def _mkelem(self, args):
#         pass
#
# class ChainExactSequence(ExactSequence, ChainComplex):
#     def __init__(self, K, subcomplex=None):
#         self.K, self.subcomplex, self.dim = K, subcomplex, K.dim
#         self._C = {}
#         ExactSequence.__init__(self, range(self.dim+2))
#         self._D = {d : Boundary(self[d], self[d-1]) for d in range(K.dim+2)}
#     def _mkelem(self, d):
#         self._C[d] = ChainGroup(self.K, d, subgroup=self._subgroup(d))
#         return self._C[d]
#
#
# class ExactSequence(Sequence):
#     def __init__(self, sequence):
#         last = self._mkzero()
#         self._sequence = []
#         self._functions = []
#         for s in sequence:
#             self._maps.append(self.function_t(last, s))
#             self._sequence.append(s)
#             last = s
#         z = self._mkzero()
#         self._maps.append(self.function_t(last, z))
#         # self._sequence.append(z)
#     def __len__(self):
#         return len(self._sequence)
#     @abstractmethod
#     def _mkzero(self):
#         pass
#     @property
#     @abstractmethod
#     def function_t(self):
#         pass
#     @property
#     @abstractmethod
#     def zero_t(self):
#         pass
#     @property
#     @abstractmethod
#     def element_t(self):
#         pass


class ChainComplex(Sequence):
    def __init__(self, it):
        self._it = it
        self._groups = {d : self._mkgroup(d) for d in self}
        self._homomorphisms = {d : self._mkhomomorphism(d) for d in self}
        assert all(self(d)(self(d+1).im).is_zero() for d in self)
    def __iter__(self):
        for d in self._it:
            yield d
    def __len__(self):
        return len(self._groups)
    def __getitem__(self, d):
        if d in self._groups:
            return self._groups[d]
        return self._mkgroup(d)
    def __call__(self, d):
        if d in self._homomorphisms:
            return self._homomorphisms[d]
        return self._mkhomomorphism(d)
    @property
    @abstractmethod
    def space_t(self):
        pass
    @abstractmethod
    def _mkgroup(self, d):
        pass
    @abstractmethod
    def _mkhomomorphism(self, d):
        pass



class ExactSequence(ChainComplex):
    def __init__(self, it):
        ChainComplex.__init__(self, it)
        assert all(self(d+1).im == self(d).ker for d in self)


class ChainBoundary(Homomorphism):
    space_t = ChainGroup
    def __init__(self, X, Y):
        assert Y.dim == X.dim - 1
        Homomorphism.__init__(self, X, Y)
    def _map(self, x):
        return x.boundary()
    def __repr__(self):
        return '-bdy->'


class SimplicialChainComplex(ChainComplex):
    space_t = ChainGroup
    def __init__(self, K, subcomplex=None, bases=None):
        self.K, self.subcomplex, self.bases, self.dim = K, subcomplex, bases, K.dim
        ChainComplex.__init__(self, tuple(reversed(range(K.dim+1))))
    def _getbasis(self, d):
        return (self.bases[d] if self.bases is not None
            and d in self.bases
            else None)
    def _mkgroup(self, d):
        return ChainGroup(self.K, d, self._getbasis(d), self._subgroup(d))
    def _subgroup(self, d):
        if self.subcomplex is not None:
            return self.subcomplex[d]
        return None
    def _mkhomomorphism(self, d):
        return ChainBoundary(self[d], self[d-1])
    def boundaries(self, d):
        return self(d+1).im
    def cycles(self, d):
        return self(d).ker
    def as_pair(self, CL):
        return SimplicialChainComplex(self.K, CL)
    def subspace(self, bases):
        return SimplicialChainComplex(self.K, self.subcomplex, bases)
    def __repr__(self):
        pairs = ['%s %s' % (str(self(d+1)), str(self[d])) for d in self]
        return '0 %s %s 0' % ('\n  '.join(pairs), self(0))


class ShortExactSequence(ExactSequence):
    def __init__(self):
        ExactSequence.__init__(self, (2,1,0))


class ChainInclusion(Monomorphism, Identity):
    space_t = ChainGroup
    def __init__(self, X, Y):
        Monomorphism.__init__(self, X, Y)
        assert all(x in self.Y for x in self.X)


class ZeroChainInclusion(ChainInclusion):
    def __init__(self, Y):
        ChainInclusion.__init__(self,  Y.zero_group(), Y)
    def __repr__(self):
        return '----->'


class ChainProjection(Epimorphism, Identity):
    space_t = ChainGroup
    def __init__(self, X, Y):
        Epimorphism.__init__(self, X, Y)
        assert all(x in self.Y or x == self.Y.zero for x in self.X)


class ZeroChainProjection(ChainProjection):
    def __init__(self, X):
        ChainProjection.__init__(self, X, X.zero_group())
    def __repr__(self):
        return '----->'

class ChainPairSES(ShortExactSequence):
    group_t = ChainGroup
    def __init__(self, G, H):
        self.G, self.H = G, H
        ShortExactSequence.__init__(self)
    def _mkhomomorphism(self, d):
        if d == 3:
            return ZeroChainInclusion(self[2])
        if d == 2:
            return ChainInclusion(self[2], self[1])
        elif d == 1:
            return ChainProjection(self[1], self[0])
        elif d == 0:
            return ZeroChainProjection(self[0])
    def _mkgroup(self, d):
        if d == 2:
            return self.H
        elif d == 1:
            return self.G
        elif d == 0:
            return self.G.as_pair(self.H)
    def __repr__(self):
        pairs = ['%s %s' % (str(self(d+1)), str(self[d])) for d in self]
        return '0 %s %s 0' % ('\n  '.join(pairs), self(0))

# ZeroChainInclusion(CL[2])
# ChainInclusion(CL[2], CK[2])
# ChainProjection(CK[2], C[2])
# ZeroChainProjection(C[2])
#
# self[dim, 2] = CL[dim]
# self[dim, 1] = CK[dim]
# self[dim, 0] = C[dim]
# self[dim] = ChainPairSES(CK[dim], CL[dim])
# self[:, 2] = SimplicialChainComplex(L)
# self[:, 1] = SimplicialChainComplex(K)
# self[:, 0] = SimplicialChainComplex(K, self[:, 2])
#
#
# self(:,2) = 0 -> SimplicialChainComplex(L)
# self(:,2) = SimplicialChainComplex(L) -> SimplicialChainComplex(K)
# self(:,1) = SimplicialChainComplex(K) -> SimplicialChainComplex(K, self[:, 2])
# self(:,0) = SimplicialChainComplex(K, self[:, 2]) -> 0
# self(dim) = ChainPairSES(CK[dim], CL[dim]) -> ChainPairSES(CK[dim-1], CL[dim-1])
# self(dim, 3) = ZeroChainInclusion(CL[dim])
# self(dim, 2) = ChainInclusion(CL[dim], CK[dim])
# self(dim, 1) = ChainInclusion(CK[dim], C[dim])
# self(dim, 0) = ZeroChainProjection(C[dim])



class ChainComplexFunctor:
    def __init__(self, X, Y):
        self.X, self.Y = X, Y
        self.im = self.Y.subspace([X(d).im for d in X])
        self.ker = self.X.subspace([X(d).ker for d in X])


class ChainComplexPairSES(ShortExactSequence):
    space_t = ChainGroup
    def __init__(self, CK, CL):
        self.CK, self.CL = CK, CL
        ShortExactSequence.__init__(self)
        self._sequences = {d : ChainPairSES(self._groups[d][2],
                                            self._groups[d][1]) for d in self}
        self._seqhomomorphisms = {d : ChainComplexFunctor(
                                        self._sequences[d],
                                        self._sequences[d-1]) for d in self}
    def _mksequence(self, dim):
        return ChainPairSES(self[dim][1], self[dim][2])
    def _mkgroup(self, idx):
        if idx == 2:
            return self.CL
        elif idx == 1:
            return self.CK
        elif idx == 0:
            return self.CK.as_pair(self.CL)
    def _mkhomomorphism(self, idx):
        return ChainComplexFunctor(self._groups[idx], self._groups[idx-1])
    def __getitem__(self, dim=None, idx=None):
        if dim is None and idx in self._groups:
            return self._groups[idx]
        elif dim in self._sequences and idx is None:
            return self._sequences[dim]
        elif dim in self._sequences and idx in self._groups:
            return self._sequences[dim][idx]
    def __call__(self, dim=None, idx=None):
        if dim is None and idx in self._groups:
            return self._groups[idx]
        elif dim in self._sequences and idx is None:
            return self._sequences[dim]
        elif dim in self._sequences and idx in self._groups:
            return self._sequences[dim][idx]
    # def _mkgroup(self, d):
    #     if d == 2:
    #         return self.CL
    #     elif d == 1:
    #         return self.CK
    #     elif d == 0:
    #         return self.CK.as_pair(self.CL)
    # def _mkhomomorphism(self, d):
    #     if d == 3:
    #         return ZeroChainInclusion(self[2])
    #     if d == 2:
    #         return ChainInclusion(self[2], self[1])
    #     elif d == 1:
    #         return ChainProjection(self[1], self[0])
    #     elif d == 0:
    #         return ZeroChainProjection(self[0])


# class ChainComplexInclusion(Monomorphism, Identity):
#     space_t = SimplicialChainComplex
#     def __init__(self, X, Y):
#         Monomorphism.__init__(self, X, Y)
#         assert all(x in self.Y for x in self.X)
#
#
# # class ZeroChainInclusion(ChainInclusion):
# #     def __init__(self, Y):
# #         ChainInclusion.__init__(self,  Y.zero_group(), Y)
# #     def __repr__(self):
# #         return '----->'
#
#
# class ChainComplexProjection(Epimorphism, Identity):
#     space_t = SimplicialChainComplex
#     def __init__(self, X, Y):
#         Epimorphism.__init__(self, X, Y)
#         assert all(x in self.Y or x == self.Y.zero for x in self.X)
#
#
# # class ZeroChainProjection(ChainProjection):
# #     def __init__(self, X):
# #         ChainProjection.__init__(self, X, X.zero_group())
# #     def __repr__(self):
# #         return '----->'
#
# class ChainComplexPairSES(ShortExactSequence):
#     def __init__(self, CK, CL):
#         self.CK, self.CL = CK, CL
#     def _mkgroup(self, d):
#         if d == 2:
#             return self.CL
#         elif d == 1:
#             return self.CK
#         elif d == 0:
#             return self.CK.as_pair(self.CL)
#     def _mkhomomorphism(self, d):
#         if d == 3:
#             return ZeroChainInclusion(self[2])
#         if d == 2:
#             return ChainInclusion(self[2], self[1])
#         elif d == 1:
#             return ChainProjection(self[1], self[0])
#         elif d == 0:
#             return ZeroChainProjection(self[0])



# class ExactSequenceZeroElement(ExactSequenceElement):
#     def __init__(self, it=None):
#         ExactSequenceElement.__init__(self, it)


# class ExactSequenceFunction(Function):
#     def __init__(self, X, Y):
#         assert (isinstance(X, ExactSequenceElement)
#             and isinstance(Y, ExactSequenceElement))
#         Function.__init__(self, X, Y)
#         assert self.im == next(self).ker
#     def __next__(self):
#         return self.__class__(self.Y, next(self.Y))
#
# class ExactSequence:
#     def __init__(self, sequence):
#         self._sequence = ExactSequenceElement(iter(sequence))
#     def __iter__(self):
#         for e in self._sequence:
#             yield e
