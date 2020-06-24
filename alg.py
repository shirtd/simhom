from collections.abc import Hashable, Set
from itertools import combinations
from sympy import Matrix, zeros
from abc import abstractmethod


from itertools import combinations
from collections.abc import Callable
from abc import abstractmethod

from simhom.util import stuple, lmap, lfilt, cref, hstack
from collections.abc import Mapping, MutableMapping
from sympy import Matrix, zeros

def csum(*C):
    return sum(*C, Chain())


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

class Complex:
    pass

class ComplexElement:
    pass

class SimplicialSkeleton(ComplexElement, MutableMapping):
    def __init__(self, dim=-1):
        self.dim, self.imap, self.S = dim, {}, []
    def add(self, s):
        assert s.dim == self.dim
        if not s in self:
            self.imap[s] = len(self)
            self.S.append(s)
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
    def __contains__(self, c):
        if isinstance(c, Simplex):
            return c in self.imap
        return all(s in self for s in c)
    def __repr__(self):
        return str(self.S)


class SimplicialComplex(Complex):
    def __init__(self, S):
        self.dim = max(len(s) for s in S) - 1
        self.S = {d : SimplicialSkeleton(d) for d in range(self.dim+1)}
        # self.imap = {d : {} for d in range(self.dim+1)}
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


class GroupElement(Hashable):
    def __init__(self, id):
        self.id = id
    def is_zero(self):
        return isinstance(self, ZeroElement)
    def __hash__(self):
        return hash(self.id)
    @abstractmethod
    def __eq__(self, e):
        pass
    @abstractmethod
    def __add__(self, g):
        pass

class ZeroElement(GroupElement):
    def __init__(self):
        GroupElement.__init__(self, tuple())
    def __add__(self, g):
        return g

class Group:
    def __init__(self, zero):
        self.zero = zero
    @abstractmethod
    def __contains__(self, g):
        pass


class FGGroup(Group):
    def __init__(self, basis, zero):
        Group.__init__(self, zero)
        self.basis = basis
    def __len__(self):
        return len(self.basis)
    def __iter__(self):
        for g in self.basis:
            yield g
    def __contains__(self, g):
        if not len(self.basis):
            return False
        M = hstack(self.to_mat(), self.to_vec(g))
        return all(p < len(self.basis) for p in M.rref()[1])
    def to_mat(self):
        return hstack(map(self.to_vec, self.basis))
    @abstractmethod
    def to_vec(self, g):
        pass
    @abstractmethod
    def to_element(self, v):
        pass


class Chain(Mapping, GroupElement):
    def __init__(self, so_pairs={}, dim=-1):
        self.omap = {s : o for s,o in so_pairs.items() if o != 0}
        self.simplices = list(sorted(self.omap))
        self.dim = self.simplices[0].dim if len(self.simplices) else dim
        assert all(s.dim == self.dim for s in self.simplices)
        id = tuple((s, self[s]*self[self[0]]) for s in self)
        GroupElement.__init__(self, id)
    def __eq__(self, c):
        return c in self and self in c
    def __add__(self, c):
        omap = {s : self[s] + c[s] for s in self.union(c)}
        if all(v == 0 for v in omap.values()):
            return ZeroChain(dim=self.dim)
        return Chain(omap, self.dim)
    def __sub__(self, c):
        return self + c * -1
    def __len__(self):
        return len(self.simplices)
    def __iter__(self):
        for s in self.simplices:
            yield s
    def __getitem__(self, s):
        if isinstance(s, Simplex):
            return self.omap[s] if s in self else 0
        return self.simplices[s]
    def __mul__(self, a):
        return Chain({s : o * a for s,o in self.items()})
    def __contains__(self, s):
        if isinstance(s, Simplex):
            return s in self.omap
        return all(ss in self and s[ss] == self[ss] for ss in s)\
            or all(ss in self and s[ss] == self[ss]*-1 for ss in s)
    def boundary(self):
        return csum(s.boundary() * o for s, o in self.items())
    def union(self, c):
        return set(self).union(set(c))
    def __lt__(self, c):
        return self.id < c.id
    def as_vec(self, imap):
        v = zeros(len(imap),1)
        for s,o in self.omap.items():
            v[imap[s]] = o
        return v
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
    def __init__(self, dim=-1, subgroup=None):
        dim = subgroup.dim if subgroup is not None else dim
        Chain.__init__(self, {}, dim)
        ZeroElement.__init__(self)
        self.basis = subgroup.basis if subgroup is not None else []
        self.subgroup = subgroup
    def __eq__(self, g):
        if g.is_zero():
            return True
        if self.subgroup is not None:
            return g in self.subgroup
    def __repr__(self):
        return '0'


class ChainGroup(FGGroup):
    def __init__(self, K, dim, basis=None, subgroup=None):
        self.K, self.dim = K, dim
        self.S = self.K[self.dim]
        if basis is None:
            basis = [Chain({s : 1}) for s in self.S]
        zero = ZeroChain(self.dim, subgroup)
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
        return ChainGroup(self.K, self.dim, basis, self.zero)
    # def __truediv__(self, H):
    #     assert H.dim == self.dim
    #     return HomologyGroup(self, H)

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
        # self.coker = self.Y.subgroup(f=lambda y: not y in self.im)
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
    def _kernel(self):
        pass


class Homomorphism(Function):
    def __init__(self, X, Y):
        Function.__init__(self, X, Y)
        assert (isinstance(X, FinitelyGeneratedGroup)
            and isinstance(Y, FinitelyGeneratedGroup)
            and is_homomorphism(self))
    def _kernel(self):
        if self.im.to_mat().is_zero:
            return self.X
        null = self.im.to_mat().nullspace()
        return self.X.subgroup(lmap(self.X.to_element, null))


class Boundary(Homomorphism):
    def __init__(self, X, Y):
        assert (isinstance(X, ChainGroup)
            and isinstance(Y, ChainGroup)
            and Y.dim == X.dim - 1)
        Homomorphism.__init__(self, X, Y)
        self.ker = self._kernel()
    def __call__(self, x):
        return x.boundary()
    def __repr__(self):
        return self.im.to_mat().__repr__()


if __name__ == '__main__':
    K = SimplicialComplex([Simplex(0,1,2),Simplex(0,2,3)])
    L = SimplicialComplex(lmap(Simplex, [[0,1],[1,2],[2,3],[0,3]]))
    CL = {d : ChainGroup(L, d) for d in range(3)}
    # CK = {d : ChainGroup(K, d, subgroup=CL[d]) for d in range(3)}
    CK = {d : ChainGroup(K, d) for d in range(3)}

    X, Y = CK[2], CK[1]

    im = ChainGroup(K, 1, [c.boundary() for c in X.basis])
    Z = hstack(im.to_vec(c) for c in Y.zero.basis)
    M = hstack(im.to_mat(), Z)
    ker = ChainGroup(K, 2, [X.to_element(m[:len(im),0]) for m in M.nullspace()])
    # M.rref()[0]
