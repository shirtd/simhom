# from simhom.simplicial import Chain
# from simhom.homology import ChainGroup
from simhom.algebra import FGGroup
from itertools import combinations
from collections.abc import Callable
from abc import abstractmethod
from simhom.util import lmap, hstack

def is_homomorphism(f):
    return all(f(x+y) == f(x)+f(y) for x,y in combinations(f.X,2))

def is_injective(f):
    return all(len(f.preimage(y)) < 2 for y in f.Y)

def is_surjective(f):
    return all(len(f.preimage(y)) > 0 for y in f.Y)

def is_bijective(f):
    return is_injective(f) and is_surjective(f)

def is_isomorphism(f):
    return is_bijective(f) and is_homomorphism(f)


# class Functor(Callable):
#     def __init__(self, A, B):
#         assert (isinstance(A, self.complex_t)
#             and isinstance(B, self.complex_t))
#         self.A, self.B = A, B
#     def __call__(self, x):
#         if isinstance(x, self.complex_t.space_t):
#             return self._map(x)
#         elif isinstance(x, self.complex_t):
#             return self._image(x)
#     # def preimage(self, y):
#     #     return {x for x in self.X if self(x) == y}
#     @property
#     @abstractmethod
#     def complex_t(self):
#         pass
#     @abstractmethod
#     def _map(self, x):
#         pass
#     @abstractmethod
#     def _preimage(self, X=None):
#         pass
#     @abstractmethod
#     def _image(self, X=None):
#         pass
#     @abstractmethod
#     def _kernel(self):
#         pass

class Function(Callable):
    def __init__(self, X, Y):
        assert (isinstance(X, self.space_t)
            and isinstance(Y, self.space_t))
        self.X, self.Y = X, Y
    def __call__(self, x):
        if isinstance(x, self.space_t.element_t):
            return self._map(x)
        elif isinstance(x, self.space_t):
            return self._image(x)
    # def preimage(self, y):
    #     return {x for x in self.X if self(x) == y}
    @property
    @abstractmethod
    def space_t(self):
        pass
    @abstractmethod
    def _map(self, x):
        pass
    @abstractmethod
    def _preimage(self, X=None):
        pass
    @abstractmethod
    def _image(self, X=None):
        pass
    @abstractmethod
    def _kernel(self):
        pass

# class Functor(Callable):
#     def __init__(self, A, B):
#
#     def __call__(self, x):
#         if isinstance(x, self.function_t):
#             return self._fun(x)
#     @abstractmethod
#     def _fun(self):
#         pass
#     @property
#     @abstractmethod
#     def function_t(self):
#         pass

# class ComplexFunctor(Monomorphism, Identity):
#     space_t = ChainGroup
#     def __init__(self, X, Y):
#         Monomorphism.__init__(self, X, Y)
#         assert all(x in self.Y for x in self.X)


class Identity(Function):
    def _map(self, x):
        return x if x in self.Y else self.Y.zero


class Homomorphism(Function):
    def __init__(self, X, Y):
        Function.__init__(self, X, Y)
        assert (isinstance(X, FGGroup)
            and isinstance(Y, FGGroup))
        for x, y in combinations(X, 2):
            try:
                assert self(x + y) == self(x) + self(y)
            except AssertionError:
                raise AssertionError('f(%s + %s) = %s != f(%s) + f(%s) = %s' %\
                                    (str(x), str(y), str(self(x + y)), str(x),
                                        str(y), str(self(x) + self(y))))
        self.im = self._image()
        self.ker = self._kernel()
        self.coker = self._cokernel()
    def _preimage(self, Y=None):
        Y = self.Y if Y is None else Y
        if isinstance(Y, self.Y.element_t):
            pre = [x for x in self.X if self(x) == Y]
        elif isinstance(Y, self.space_t):
            pre = [x for x in self.X if self(x) in Y]
        else:
            raise TypeError('%s is not compatable with types %s, %s' %\
                    (str(Y), str(self.space_t), str(self.Y.element_t)))
        return self.X.subgroup(pre)
    def _image(self, X=None):
        X = self.X if X is None else X
        return self.Y.subgroup(lmap(self, X))
    def _kernel(self):
        if self.im.is_zero():
            return self.X
        im = lmap(self, self.X) + list(self.Y.zero.basis)
        null = hstack(self.Y.to_vec(c) for c in im).nullspace()
        ker = [self.X.to_element(m[:len(self.X),0]) for m in null]
        return self.X.subgroup(ker)
    def _cokernel(self):
        coker = [y for y in self.Y if not y in self.im]
        return self.Y.subgroup(coker)

class Monomorphism(Homomorphism):
    def __init__(self, X, Y):
        Homomorphism.__init__(self, X, Y)
        # for y in self.Y:
        #     pre = self._preimage(y)
        #     try:
        #         assert len(pre) < 2
        #     except AssertionError:
        #         raise AssertionError(self._error_str(y, pre))
    def _error_str(self, y, pre):
        return '%s\n\t%s' % ('Function is not injective:', '%s = %s'\
                % (' = '.join(['f(%s)' % str(x) for x in pre]), str(y)))
    def __repr__(self):
        return '-inj->'

class Epimorphism(Homomorphism):
    def __init__(self, X, Y):
        Homomorphism.__init__(self, X, Y)
        # try:
        #     assert self.im == self.Y
        # except AssertionError:
        #     raise AssertionError(self._error_str())
    def _error_str(self):
        return '%s\n\t%s' % ('Function is not surjective:', 'f(%s) = %s != %s'\
                            % (str(self.X), str(self.im), str(self.Y)))
    def __repr__(self):
        return '-sur->'


class Isomorphism(Epimorphism, Monomorphism):
    def __init__(self, X, Y):
        Epimorphism.__init__(self, X, Y)
        Monomorphism.__init__(self, X, Y)
    def __repr__(self):
        return '-iso->'


# class Map:
#     def __init__(self, map, X, Y):
#         self._map, self.X, self.Y = map, X, Y
#     def __call__(self, x):
#         return self._map(x)
#     def preimage(self, y):
#         preim = [x for x in self.X if self(x) == y]
#         return self.X.subspace(preim)
#     def image(self):
#         im = [self(x) for x in self.X.basis]
#         # return self.Y.subspace(im)
#         return im
#     def kernel(self):
#         ker = [x for x in self.X.basis if self(x) == self.Y.zero]
#         # return self.Y.subspace(ker)
#         return ker
#     def cokernel(self):
#         coker = [y for y in self.Y.basis if not y in self.image()]
#         # return self.Y.subspace(coker)
#         return coker
#
# # class ChainMap(Map):
# #     def __init__(self, map, X, Y):
# #         assert isinstance(X, ChainGroup) and isinstance(Y, ChainGroup)
# #         Map.__init__(self, map, X, Y)
#
# class Inclusion(Map):
#     def __init__(self, X, Y):
#         map = lambda x: x
#         Map.__init__(self, map, X, Y)
#
# class InducedHomomorphism(Map):
#     def __init__(self, f, X, Y):
#         self.CX, self.CY = f.X, f.Y
#         self.chain_map = f
#         def map(l):
#             fcs = [self.chain_map(self.CX.to_chain(v)) for v in self.X[l]]
#             Fcs = {Y(fc) for fc in fcs}
#             assert len(Fcs) == 1
#             return next(iter(Fcs))
#         Map.__init__(self, map, X, Y)
