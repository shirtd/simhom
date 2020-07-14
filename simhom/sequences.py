from __future__ import annotations


from simhom.homology import *
from simhom.util import lmap, ker, im, pad
import numpy as np

# class ChainComplexIterator:
#     @abstractmethod
#     def __iter__(self):
#         pass
#     @abstractmethod
#     def __next__(self):
#         pass

class ChainComplex:
    def __init__(self, seq, mkzero, idx=0):
        self._mkzero = mkzero
        self._seq = seq
        self.idx = idx
        self._maps = {d : self[d] >> self[d-1] for d in self}
        # assert all(self(d)(self(d+1).im) == 0 for d in self)
    def __len__(self):
        return len(self._seq)
    def __getitem__(self, d):
        if d in self._seq:
            return self._seq[d]
        return self._mkzero(d)
    def __call__(self, d):
        if d in self._maps:
            return self._maps[d]
        return self[d] >> self[d-1]
    def __eq__(self, other):
        return (other == 0
            and all(self[d] == 0 for d in self))
    @property
    @abstractmethod
    def space_t(self):
        pass
    @abstractmethod
    def __iter__(self):
        pass
    @abstractmethod
    def get_idx(self, U):
        pass
    def __repr__(self):
        pairs = ['%s %s' % (str(self(d+1)), str(self[d])) for d in self]
        return 'O %s %s O' % ('\n  '.join([pad(s, 9) for s in pairs]), self(0))

class SimplicialChainComplex(ChainComplex):
    space_t = ChainGroup
    @classmethod
    def simplicial_init(cls, K, idx=0):
        seq = {d : ChainGroup.simplicial_init(K, d) for d in range(K.dim+1)}
        return SimplicialChainComplex(seq, K.dim, idx=idx)
    def __init__(self, seq, dim, mkzero=None, idx=0):
        self.dim = dim
        if mkzero is None:
            mkzero = lambda d: ChainGroup([], ZeroChain(dim=d))
        ChainComplex.__init__(self, seq, mkzero, idx)
    def get_idx(self, U):
        return U.dim
    def __iter__(self):
        for d in reversed(range(self.dim+1)):
            yield d
    def subcomplex(self, seq, mkzero=None, idx=None):
        mkzero = self._mkzero if mkzero is None else mkzero
        idx = self.idx if idx is None else idx
        return SimplicialChainComplex(seq, self.dim, mkzero, idx)
    def boundaries(self, d=None):
        if d is None:
            seq = {d : self(d+1).im for d in self}
            return self.subcomplex(seq)
        return self(d+1).im
    def cycles(self, d=None):
        if d is None:
            seq = {d : self(d).ker for d in self}
            return self.subcomplex(seq)
        return self(d).ker
    def __truediv__(self, other):
        seq = {d : self[d] / other[d] for d in self}
        mkzero = lambda d: ChainQuotient(self._mkzero(d), self._mkzero(d))
        return self.subcomplex(seq, mkzero)
    def __rshift__(self, other):
        return InducedComplexMap(self, other)


class ComplexMap(Morphism[ChainComplex]):
    def __init__(self, X: ChainComplex, Y: ChainComplex) -> None:
        Morphism.__init__(self, X, Y)
    def __call__(self, x: Union[Group, ChainComplex])\
            -> Union[Group, ChainComplex]:
        if isinstance(x, Group):
            return self._group_map(x)
        elif isinstance(x, ChainComplex):
            return self._complex_map(x)
        else:
            raise TypeError(str(x))
    @abstractmethod
    def _group_map(self, x: Group) -> Group:
        pass
    @abstractmethod
    def _complex_map(self, U: ChainComplex) -> ChainComplex:
        pass


class InducedComplexMap(ComplexMap):
    group_t, complex_t = ChainGroup, ChainComplex
    def __init__(self, X, Y):
        ComplexMap.__init__(self, X, Y)
        self._f = {d : X[d] >> Y[d] for d in X}
        self.im = self(self.X)
        self.ker = self._kernel()
    def _group_map(self, U):
        idx = self.X.get_idx(U)
        if not idx in self._f:
            return (self.X[idx] >> self.Y[idx])(U)
        return self._f[self.X.get_idx(U)](U)
    def _complex_map(self, C):
        seq = {d : self(C[d]) for d in C}
        return self.Y.subcomplex(seq)
    def _kernel(self, zero_group=None):
        seq = {d : self._f[d].ker for d in self.X}
        return self.X.subcomplex(seq)
    def __repr__(self):
        return '|&&&&>'

class ExactSequence(ChainComplex):
    def __init__(self, seq, mkzero):
        ChainComplex.__init__(self, seq, mkzero)
        # assert all(self(d+1).im == self(d).ker for d in self)


class ShortExactSequence(ExactSequence):
    def __init__(self, seq, mkzero):
        ExactSequence.__init__(self, seq, mkzero)
    def __iter__(self):
        for d in (2, 1, 0):
            yield d

# class ChainPairSES(ShortExactSequence):
#     space_t = ChainGroup
#     @classmethod
#     def simplicial_init(cls, K, L, dim):
#         CL = ChainGroup.simplicial_init(L, dim)
#         CK = ChainGroup.simplicial_init(K, dim)
#         return ChainPairSES.group_init(CK, CL, dim)
#     @classmethod
#     def group_init(cls, CK, CL, dim):
#         seq = {2 : CL, 1 : CK, 0 : CK / CL}
#         seq[2].idx, seq[1].idx, seq[0].idx = 2, 1, 0
#         return ChainPairSES(seq, dim)
#     def __init__(self, seq, dim, mkzero=None):
#         self.dim = dim
#         if mkzero is None:
#             mkzero = lambda d: ChainGroup([], ZeroChain(dim=dim))
#         ShortExactSequence.__init__(self, seq, mkzero)
#     def get_idx(self, U):
#         return U.idx
#     def subcomplex(self, seq, mkzero=None):
#         mkzero = self._mkzero if mkzero is None else mkzero
#         return ChainPairSES(seq, self.dim, mkzero)
#     def __truediv__(self, other):
#         seq = {d : self[d] / other[d] for d in self}
#         mkzero = lambda d: ChainQuotient(self._mkzero(d), self._mkzero(d))
#         return self.subcomplex(seq, mkzero)
#     def __rshift__(self, other):
#         return InducedComplexMap(self, other)

class ChainComplexPairSES(ShortExactSequence):
    space_t = ChainComplex
    @classmethod
    def simplicial_init(cls, K, L):
        CL = SimplicialChainComplex.simplicial_init(L, 2)
        CK = SimplicialChainComplex.simplicial_init(K, 1)
        return ChainComplexPairSES.group_init(CK, CL)
    @classmethod
    def group_init(cls, CK, CL):
        seq = {2 : CL, 1 : CK, 0 : CK / CL}
        return ChainComplexPairSES(seq, CK.dim)
    def __init__(self, seq, dim, mkzero=None):
        self.dim = dim
        if mkzero is None:
            mkzero = lambda d: SimplicialChainComplex({}, dim)
        ShortExactSequence.__init__(self, seq, mkzero)
    def get_idx(self, U):
        return U.idx
    def subcomplex(self, seq, mkzero=None):
        mkzero = self._mkzero if mkzero is None else mkzero
        return ChainComplexPairSES(seq, self.dim, mkzero)
    def __truediv__(self, other):
        seq = {d : self[d] / other[d] for d in self}
        mkzero = lambda d: self._mkzero(d) / other._mkzero(d)
        return self.subcomplex(seq, mkzero)
    def __rshift__(self, other):
        return InducedComplexMap(self, other)
    def cycles(self, idx=None, dim=None):
        if idx is None and dim is None:
            seq = {i : self[i].cycles() for i in self}
            # seq[2].idx, seq[1].idx, seq[0].idx = 2, 1, 0
            return ChainComplexPairSES(seq, self.dim)
        elif idx is None:
            seq = {i : self[i].cycles(dim) for i in self}
            # seq[2].idx, seq[1].idx, seq[0].idx = 2, 1, 0
            return ChainPairSES(seq, self.dim)
        elif dim is None:
            return self[idx].cycles()
        return self[idx].cycles(dim)
    def boundaries(self, idx=None, dim=None):
        if idx is None and dim is None:
            seq = {i : self[i].boundaries() for i in self}
            # seq[2].idx, seq[1].idx, seq[0].idx = 2, 1, 0
            return ChainComplexPairSES(seq, self.dim)
        elif idx is None:
            seq = {i : self[i].boundaries(dim) for i in self}
            # seq[2].idx, seq[1].idx, seq[0].idx = 2, 1, 0
            return ChainPairSES(seq, self.dim)
        elif dim is None:
            return self[idx].boundaries()
        return self[idx].boundaries(dim)




class HomologyPairLES(ExactSequence):
    @classmethod
    def complex_init(self, C):
        seq = {3*d + i : C[i][d] for d in range(C.dim+1) for i in C}
        for k in seq:
            seq[k].idx = k
        def mkzero(d):
            return ChainGroup([], ZeroChain(dim=d))\
                    / ChainGroup([], ZeroChain(dim=d))
        return HomologyPairLES(seq, C.dim, mkzero)
    def __init__(self, seq, dim, mkzero):
        self.dim = dim
        ExactSequence.__init__(self, seq, mkzero)
    def get_idx(self, U):
        return U.idx
    def __iter__(self):
        for i in reversed(range((self.dim+1)*3)):
            yield i
    def __getitem__(self, t):
        t = 3*t[0] + t[1] if isinstance(t, tuple) else t
        return ExactSequence.__getitem__(self, t)
    def __call__(self, *t):
        t = 3*t[0] + t[1] if len(t) > 1 else t[0]
        return ExactSequence.__call__(self, t)
    def subcomplex(self, seq, mkzero=None):
        mkzero = self._mkzero if mkzero is None else mkzero
        return HomologyPairLES(seq, self.dim, mkzero)
    def __rshift__(self, other):
        return InducedComplexMap(self, other)
