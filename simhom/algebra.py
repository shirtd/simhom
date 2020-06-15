from sympy import Matrix, zeros
from simhom.util import hstack, stuple
from simhom.simplicial import Chain


class Group:
    def __init__(self, basis, zero, op=lambda uv: uv[0]+uv[1]):
        self.basis, self.zero, self.op = basis, zero, op
        self.elements = self.basis + [self.zero]
    def __iter__(self):
        # for v in self.basis:
        for v in self.elements:
            yield v
    def __len__(self):
        return len(self.basis)


class VectorSpace(Group):
    def __init__(self, basis, zero=None):
        self.m = next(iter(basis)).shape[0] if len(basis) else 0
        assert all(u.shape[0] == self.m for u in basis)
        zero = zeros(self.m, 1) if zero is None else zero
        # basis = [v for v in basis if not v == zero]
        Group.__init__(self, basis, zero)
        self.hash = stuple(map(tuple, self.cref().columnspace()))
    def __truediv__(self, V):
        return Quotient(self, V)
    def to_vec(self, u):
        return u
    def to_mat(self):
        if all(u.shape[0] == 0 for u in self):
            return zeros(0, len(self))
        return hstack(*[self.to_vec(v) for v in self.basis])
    def cref(self):
        return self.to_mat().T.rref()[0].T
    def rref(self):
        return self.to_mat().rref()[0]
    def __contains__(self, u):
        if not len(self):
            return False
        M = hstack(self.to_mat(), self.to_vec(u))
        return all(p < len(self) for p in M.rref()[1])
        # return M.T.rref()[0].T[:,-1].is_zero
    def __repr__(self):
        return self.to_mat().__repr__()
    def __eq__(self, V):
        if self.to_mat().is_zero and V.to_mat().is_zero:
            return True
        return all(u in V for u in self) and all(v in self for v in V)
    def __hash__(self):
        return hash(self.hash)
    def subspace(self, V):
        assert all(v in self for v in V)
        return VectorSpace(V, self.zero)


class Coset(VectorSpace):
    def __init__(self, u, V):
        uV = hstack(*[u + v for v in V]).T.rref()[0].T
        basis = [uV[:,i] for i in range(uV.shape[1])]
        VectorSpace.__init__(self, basis, V.cref())


class Quotient(Group):
    def __init__(self, U, V):
        self.U, self.V = U, V
        zero = Coset(U.zero, V)
        self.classes = {}
        for u in U:
            c = Coset(u, V)
            if not c in self.classes:
                self.classes[c] = []
            self.classes[c].append(u)
        def op(A, B):
            reps, AB = [], set()
            for a in self.classes[A]:
                for b in self.classes[B]:
                    ab = self.U.op((a, b))
                    AB.add(Coset(ab, self.V))
                    reps.append(ab)
            assert len(AB) == 1
            return next(iter(AB))
        basis = [c for c in self.classes.keys() if not c == zero]
        Group.__init__(self, basis, zero, op)
    def __getitem__(self, l):
        return self.classes[l]
    def __call__(self, u):
        if isinstance(u, Chain):
            u = self.U.to_vec(u)
        return Coset(u, self.V)
