from simhom import *

if __name__ == '__main__':
    S = [[0,1,2]]
    K = SimplicialComplex(S)
    L = SimplicialComplex([[0,1], [1,2], [0,2]])
    CK, CL = ChainComplex(K), ChainComplex(L)
    E = LongExactSequence(CK, CL)

    f = Inclusion(CL[0], CK[0])
    F = InducedHomomorphism(f, E.HL[0], E.HK[0])
    l = E.HL[0].basis[0]
    print(F(l))
    c = next(iter(E.HL[1].classes.values()))[0]
    # c = F.CX.to_chain(F.X[list(F.X)[0]][0])
