from simhom.simplicial import Chain
from simhom.homology import Homology, ChainGroup, ChainComplex

class Map:
    def __init__(self, map, X, Y):
        self._map, self.X, self.Y = map, X, Y
    def __call__(self, x):
        return self._map(x)
    def preimage(self, y):
        preim = [x for x in self.X if self(x) == y]
        return self.X.subspace(preim)
    def image(self):
        im = [self(x) for x in self.X]
        return self.Y.subspace(im)
    def kernel(self):
        ker = [x for x in self.X if self(x) == self.Y.zero]
        return self.Y.subspace(ker)
    def cokernel(self):
        coker = [y for y in self.Y if not y in self.image()]
        return self.Y.subspace(coker)

# class ChainMap(Map):
#     def __init__(self, map, X, Y):
#         assert isinstance(X, ChainGroup) and isinstance(Y, ChainGroup)
#         Map.__init__(self, map, X, Y)

class Inclusion(Map):
    def __init__(self, X, Y):
        map = lambda x: x
        Map.__init__(self, map, X, Y)

class InducedHomomorphism(Map):
    def __init__(self, f, X, Y):
        self.CX, self.CY = f.X, f.Y
        self.chain_map = f
        def map(l):
            fcs = [self.chain_map(self.CX.to_chain(v)) for v in self.X[l]]
            Fcs = {Y(fc) for fc in fcs}
            assert len(Fcs) == 1
            return next(iter(Fcs))
        Map.__init__(self, map, X, Y)
