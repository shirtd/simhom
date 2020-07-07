from sympy import Matrix, zeros

def stuple(s, *args, **kw):
    return tuple(sorted(s, *args, **kw))

def slist(s, *args, **kw):
    return list(sorted(s, *args, **kw))

def hstack(*L):
    try:
        return Matrix([l.T for l in L]).T
    except AttributeError:
        return Matrix([l.T for l in iter(*L)]).T

def vstack(*L):
    return Matrix(L)

def pad(s, pad=0):
    return s.replace('\n', '\n%s' % (' '*pad))

# def rref(M): return M.rref()[0]
#
# def cref(M): return rref(M.T).T

def tmap(f, t): return tuple(map(f, t))
def stmap(f, t): return stuple(map(f, t))
def lmap(f, l): return list(map(f, l))
def smap(f, s): return set(map(f, s))
def tfilt(f, t): return tuple(filter(f, t))
def stfilt(f, t): return stuple(filter(f, t))
def lfilt(f, l): return list(filter(f, l))
def sfilt(f, s): return set(filter(f, s))

def to_mat(basis):
    elems = list({e for b in basis for e in b})
    imap = {s : i for i, s in enumerate(elems)}
    vs = []
    for c in basis:
        v = zeros(len(elems),1)
        for s in c:
            v[imap[s]] = int(c[s])
        vs.append(v)
    return hstack(vs)


def ker(f): return f.ker
def im(f): return f.im
