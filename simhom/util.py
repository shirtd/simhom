from sympy import *

def stuple(s):
    return tuple(sorted(s))

def csum(*C):
    return sum(*C, Chain())

def hstack(*L):
    return Matrix([l.T for l in L]).T

def vstack(*L):
    return Matrix(L)

def pad(s, pad=0):
    return s.replace('\n', '\n%s' % (' '*pad))

def lmap(f, v):
    return list(map(f, v))
