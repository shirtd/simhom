from simhom.algebra import BasisElement, ZeroElement
from simhom.util import stuple, lmap, lfilt, cref, slist
from collections.abc import Mapping, MutableMapping
from sympy import Matrix, zeros
from abc import abstractmethod

def csum(*C):
    return sum(*C, ZeroChain())


''''''''''''''''''
''' PRIMITIVES '''
''''''''''''''''''


class Simplex(tuple):
    def __new__(cls, *v):
        try:
            return tuple.__new__(cls, stuple(*iter(v)))
        except:
            return tuple.__new__(cls, stuple(v))
    def __init__(self, *v):
        tuple.__init__(self)
        self.dim = len(self) - 1
    def boundary(self):
        if self.dim == 0:
            return Chain()
        return Chain({self.get_face(i) : (-1)**i for i in range(self.dim+1)})
    def get_face(self, i):
        return Simplex(self[:i]+self[i+1:])
    def __contains__(self, s):
        if isinstance(s, Simplex) or isinstance(s, Chain):
            return all(v in self for v in s)
        return s in self


class Chain(BasisElement):
    atom_t = Simplex
    def __init__(self, so_pairs={}, dim=-1):
        map = {s : o for s, o in so_pairs.items() if o != 0}
        BasisElement.__init__(self, map)
        self.id = tuple((s, o*self.id[0][1]) for s,o in self.id)
        self.dim = self.atoms[0].dim if len(self) else dim
        try:
            assert all(s.dim == self.dim for s in self)
        except AssertionError as e:
            print(self)
            raise e
    def _mkchain(self, map, dim=-1):
        if all(v == 0 for v in map.values()) or not len(map):
            return ZeroChain(dim=dim)
        return Chain({s : o for s,o in map.items()}, dim)
    def __eq__(self, h):
        return (BasisElement.__eq__(self, h)
            or BasisElement.__eq__(self, ~h))
    def __hash__(self):
        return BasisElement.__hash__(self)
    def __getitem__(self, s):
        if isinstance(s, Simplex):
            return BasisElement.__getitem__(self, s)
        return self.simplices[s]
    def __add__(self, c):
        simplices = set(self).union(set(c))
        map = {s : self[s] + c[s] for s in simplices}
        return self._mkchain(map, self.dim)
    def __mul__(self, a):
        map = {s : o * a for s,o in self.items()}
        return self._mkchain(map, self.dim)
    def __invert__(self):
        return self * -1
    def boundary(self):
        return csum(s.boundary() * o for s, o in self.items())
    def __repr__(self):
        s = ''
        for i, (t, o) in enumerate(self.items()):
            if i > 0 and o > 0:
                s += ' + '
            if o < 0:
                s += ' - '
            if abs(o) > 1:
                s += '%d*' % abs(o)
            s += str(t)
        return s


class ZeroChain(ZeroElement, Chain):
    def __init__(self, subgroup=None, dim=-1):
        dim = subgroup.dim if subgroup is not None else dim
        Chain.__init__(self, {}, dim)
        ZeroElement.__init__(self, subgroup)
    def __repr__(self):
        return '0'


class ComplexElement:
    # def __init__(self, left, right):
    #     self._left, self._right = left, right
    #     self._from = morphism_t(left, self)
    #     self._to = morphism_t(self, right)
    # # def __next__(self):
    # #     return self._right
    # @abstractmethod
    # def __rshift__(self, sup):
    #     pass
    @property
    @abstractmethod
    def morphism_t(self):
        pass
    @property
    @abstractmethod
    def element_t(self):
        pass

# class ComplexMorphism(Function):
#     def __init__(self, X, Y):
#         Function.__init__(self, X, Y)
#         self.im = self.Y>>self
#         assert not len(next(self.im))
#     def __iter__(self):
#         for x in self.X:
#             yield self(x)
#
# class SimplicialBoundary(ComplexMorphism):
#     def __init__(self, X, Y):
#         assert (isinstance(X, SimplicialSkeleton)
#             and isinstance(Y, SimplicialSkeleton)
#             and Y.dim == X.dim - 1)
#         ComplexMorphism.__init__(self, X, Y)
#     def __call__(self, x):
#         return {s for s in x.boundary()}

class SimplicialSkeleton:#(ComplexElement, MutableMapping):
    # morphism_t = Boundary
    def __init__(self, dim=-1):#, S=[]):
        self.dim, self.S, self.imap = dim, [], {}
        # for s in S:
        #     self.add(s)
    # def __rshift__(self, sup):
    #     S = {s for s in S for S in sup}
    #     return SimplicialSkeleton(S=slist(S, key=lambda s: self[s]))
    def __len__(self):
        return len(self.S)
    def __getitem__(self, s):
        if isinstance(s, Simplex):
            return self.imap[s]
        return self.S[s]
    def __setitem__(self, s, i):
        self.imap[s] = i
    def __delitem__(self, s):
        del self.imap[s]
        self.S.remove(s)
    def __iter__(self):
        for s in self.S:
            yield s
    def to_vec(self, c):
        assert c.dim == self.dim
        v = zeros(len(self),1)
        if isinstance(c, Simplex):
            v[self[c]] = 1
        else:
            for s, o in c.items():
                v[self[s]] = o
        return v
    def add(self, s):
        if self.dim == -1:
            self.dim = s.dim
        elif self.dim != s.dim:
            raise Exception('dimension of %s does not match' % str(s))
        if not s in self:
            self.imap[s] = len(self)
            self.S.append(s)
    def __contains__(self, c):
        if isinstance(c, Simplex):
            return c in self.imap
        return all(s in self for s in c)
    def __repr__(self):
        return str(self.S)

''''''''''''''''''
''' STRUCTURES '''
''''''''''''''''''


class Complex:
    # def __init__(self, sequence):
    #     assert all(next(next(s)).is_zero() for s in sequence) == Zero
    @property
    @abstractmethod
    def element_t(self):
        pass


class SimplicialComplex:#(Complex):
    # element_t = SimplicialSkeleton
    def __init__(self, S):
        self.dim = max(len(s) for s in S) - 1
        self.S = {d : SimplicialSkeleton(d) for d in range(self.dim+1)}
        for s in S:
            self.add(Simplex(s))
    def to_vec(self, c):
        if c.dim < 0:
            return zeros(0,1)
        return self[c.dim].to_vec(c)
    def add(self, s):
        self[s.dim].add(s)
        if s.dim > 0:
            for t in s.boundary():
                self.add(t)
    def get_boundaries(self, dim):
        return [self.to_vec(s.boundary()) for s in self[dim]]
    def __getitem__(self, d):
        if d > self.dim or d < 0:
            return SimplicialSkeleton(d)
        return self.S[d]
    def __iter__(self):
        for k,v in self.S.items():
            yield (k,v)
    def __contains__(self, c):
        return c in self[c.dim]
    def __repr__(self):
        return '{%s}' % ',\n '.join(['%d: %s' % (d, str(s)) for d,s in self])
