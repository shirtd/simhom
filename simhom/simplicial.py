from itertools import combinations
from collections.abc import Iterable
from scipy.linalg import lu
from sympy import *
import numpy as np
import numpy

from simhom.util import stuple, csum

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
        return Chain((self.get_face(i), (-1)**i) for i in range(self.dim+1))
    def get_face(self, i):
        return Simplex(self[:i]+self[i+1:])
    def __contains__(self, s):
        if isinstance(s, Simplex) or isinstance(s, Chain):
            return all(v in self for v in s)
        return s in self


class Chain:
    def __init__(self, so_pairs=[]):
        self.omap = {s : o for s,o in so_pairs if o != 0}
        self.simplices = list(sorted(self.omap.keys()))
        if len(self):
            self.dim = next(iter(self)).dim
            assert all(s.dim == self.dim for s in self)
        else:
            self.dim = -1
        self.is_zero = len(self) == 0
    def boundary(self):
        return csum(s.boundary() * o for s, o in self.items())
    def __len__(self):
        return len(self.simplices)
    def orient(self, s):
        return self.omap[s] if s in self.omap else 0
    def items(self):
        for k,v in self.omap.items():
            yield (k,v)
    def __getitem__(self, s):
        if isinstance(s, Simplex):
            return self.omap[s] if s in self else 0
        return self.simplices[s]
    def __iter__(self):
        for s in self.simplices:
            yield s
    def union(self, c):
        return set(self).union(set(c))
    def __add__(self, c):
        return Chain((s, self[s] + c[s]) for s in self.union(c))
    def __mult__(self, a):
        return Chain((s, o * a) for s,a in self.items())
    def __contains__(self, s):
        if isinstance(s, Simplex):
            return s in self.omap
        return all(ss in self and s[ss] == self[ss])\
            or all(ss in -self and s[ss] == self[ss])
    def __eq__(self, c):
        return c in self and self in c
    def __repr__(self):
        if self.is_zero:
            return '0'
        s = ''
        for i, (t, o) in enumerate(self.items()):
            if i > 0 and o > 0:
                s += ' + '
            if o < 0:
                s += ' - '
            s += str(t)
        return s


class SimplicialComplex:
    def __init__(self, S):
        self.dim = max(len(s) for s in S) - 1
        self.S = {d : list() for d in range(self.dim+1)}
        self.imap = {d : {} for d in range(self.dim+1)}
        for s in S:
            self.add(Simplex(s))
    def to_vec(self, c):
        if isinstance(c, Simplex):
            v = zeros(len(self[c.dim]),1)
            v[self.imap[c.dim][c]] = 1
            return v
        if c.dim < 0:
            return zeros(0,1)
        v = zeros(len(self[c.dim]), 1)
        for s, o in c.items():
            v[self.imap[s.dim][s]] = o
        return v
    def closure(self, s):
        if s.dim > 0:
            for t in s.boundary():
                self.add(t)
    def add(self, s):
        if not s in self[s.dim]:
            self.imap[s.dim][s] = len(self.S[s.dim])
            self.S[s.dim].append(s)
        self.closure(s)
    def get_boundaries(self, dim):
        return [self.to_vec(s.boundary()) for s in self[dim]]
    def __getitem__(self, d):
        if d > self.dim or d < 0:
            return list()
        return self.S[d]
    def __iter__(self):
        for k,v in self.S.items():
            yield (k,v)
    def __repr__(self):
        return '{%s}' % ',\n '.join(['%d: %s' % (d, str(s)) for d,s in self])
    def __contains__(self, c):
        return c in self[c.dim]
