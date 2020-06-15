from simhom.algebra import Group, VectorSpace, Coset, Quotient
from simhom.simplicial import Chain
from simhom.util import pad
from sympy import eye, zeros

class ChainGroup(VectorSpace):
    def __init__(self, K, dim, basis=None, zero=None):
        self.K, self.dim = K, dim
        basis = eye(len(K[dim])).columnspace() if basis is None else basis
        if not len(basis) and zero is None:
            zero = zeros(len(K[dim]),1)
        VectorSpace.__init__(self, basis, zero)
    def to_vec(self, c):
        if isinstance(c, Chain):
            return self.K.to_vec(c)
        return super().to_vec(c)
    def as_subspace(self, C):
        basis = [C.to_vec(self.to_chain(v)) for v in self.basis]
        return ChainGroup(C.K, self.dim, basis)
    def subspace(self, C):
        basis = [self.to_vec(C.to_chain(v)) for v in C.basis]
        return ChainGroup(self.K, self.dim, basis, self.zero)
    def restrict(self, chains):
        return
    def to_chain(self, v):
        return Chain((self.K[self.dim][i], e) for i,e in enumerate(v) if e != 0)
    def __repr__(self):
        return 'span(%s)' % ', '.join([str(self.to_chain(v)) for v in self.basis])
    def __contains__(self, c):
        if isinstance(c, Chain):
            c = self.to_vec(c)
        return super().__contains__(c)


class Boundary:
    def __init__(self, K, dim):
        self.K, self.dim = K, dim
        self.im = ChainGroup(K, dim-1, K.get_boundaries(dim))
        self.ker = ChainGroup(K, dim, self.im.to_mat().nullspace())
    def to_vec(self, c):
        return self.K.to_vec(c)
    def __repr__(self):
        return self.im.__repr__()


class RelativeBoundary(Boundary):
    def __init__(self, CK, CL, dim):
        self.K, self.L, self.dim = CK.K, CL.K, dim
        CL_K = CL[dim-1].as_subspace(CK[dim-1])
        self.im = ChainGroup(self.K, dim-1, CK.D[dim].im.basis + CL_K.basis)
        sup = [self.K.to_vec(s) for s in self.K[dim] if s.boundary() in CL_K]
        self.ker = ChainGroup(self.K, dim, CK.D[dim].ker.basis + sup)


class ChainComplex:
    def __init__(self, K):
        self.K, self.dim = K, K.dim
        self.C = {d : ChainGroup(K, d) for d,_ in K}
        self.D = {d : Boundary(K, d) for d,_ in K}
    def cycles(self, dim):
        if dim in self.D:
            return self.D[dim].ker
        return Boundary(self.K, dim).ker
    def boundaries(self, dim):
        if dim+1 in self.D:
            return self.D[dim+1].im
        return Boundary(self.K, dim+1).im
    def __getitem__(self, d):
        if d in self.C:
            return self.C[d]
        return ChainGroup(self.K, d)


class RelativeChainComplex(ChainComplex):
    def __init__(self, CK, CL):
        self.C, self.CK, self.CL = CK.C, CK, CL
        self.K, self.L, self.dim = CK.K, CL.K, CK.dim
        self.D = {d : RelativeBoundary(CK, CL, d) for d in range(self.dim+1)}


class Homology(Quotient):
    def __init__(self, C, dim):
        self.C, self.dim = C, dim
        Quotient.__init__(self, C.cycles(dim), C.boundaries(dim))
        # self.classes = {l : lmap(C[dim].to_chain, v) for l,v in self.classes.items()}
        # self.cmap = {}
        # for l, c in self.classes.items():
        #     self.cmap[self.C[dim].to_chain(c)] = l
    def print_class(self, l):
        cs = self.classes[l]
        crep = [self.C[self.dim].to_chain(c) for c in cs if not c.is_zero]
        s = ['[%s]' % str(c) for c in crep]
        s = ['[0]'] + s if l == self.zero else s
        # s = ['[%s]' % str(c) for c in self.classes[l]]
        ss = '%s%s' % (s[0], (' = %s' % s[1]) if len(s) > 1 else '')
        return ('\n%s = ' % (' ' * len(s[0]))).join([ss] + s[2:])
    # def get_class(self, c):
    #     return self.cmap[c]
    def __repr__(self):
        cstr = [self.print_class(l) for l in self]
        return 'FAB(%s)' % pad(',\n'.join(cstr), 4)


class LongExactSequence:
    def __init__(self, CK, CL):
        self.CK, self.CL = CK, CL
        self.K, self.L, self.dim = CK.K, CL.K, CK.dim
        self.C = RelativeChainComplex(CK, CL)
        self.HK = {d : Homology(self.CK, d) for d in range(CK.dim+1)}
        self.HL = {d : Homology(self.CL, d) for d in range(CL.dim+1)}
        self.HKL = {d : Homology(self.C, d) for d in range(self.dim+1)}
