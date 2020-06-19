# from simhom.simplicial import Chain
# from simhom.homology import ChainGroup
from simhom.algebra import FinitelyGeneratedGroup
from itertools import combinations
from collections.abc import Callable
from abc import abstractmethod
from simhom.util import lmap

def is_homomorphism(f):
    return all(f(x+y) == f(x)+f(y) for x,y in combinations(f.X,2))

def is_injective(f):
    return all(len(f.preimage(y)) < 2 for y in self.Y)

def is_surjective(f):
    return all(len(f.preimage(y) > 0 for y in self.Y))

def is_bijective(f):
    return is_injective(f) and is_surjective(f)

def is_isomorphism(f):
    return is_bijective(f) and is_homomorphism(f)


class Function(Callable):
    def __init__(self, X, Y):
        self.X, self.Y = X, Y
        self.im = self.image(self.X)
        self.coker = self.Y.subgroup(f=lambda y: not y in self.im)
    def image(self, x=None):
        x = self.X if x is None else x
        if self.X.is_element(x):
            return self.image({x})
        return self.Y.subgroup(lmap(self, x))
    def preimage(self, y=None):
        y = self.Y if y is None else y
        if self.Y.is_element(y):
            return self.preimage({y})
        return self.X.subgroup(f=lambda x: self(x) in y)
    @abstractmethod
    def kernel(self):
        pass


class Homomorphism(Function):
    def __init__(self, X, Y):
        Function.__init__(self, X, Y)
        assert (isinstance(X, FinitelyGeneratedGroup)
            and isinstance(Y, FinitelyGeneratedGroup)
            and is_homomorphism(self))
    def kernel(self):
        if self.im.to_mat().is_zero:
            return self.X
        null = self.im.to_mat().nullspace()
        return self.X.subgroup(lmap(self.X.to_element, null))


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
