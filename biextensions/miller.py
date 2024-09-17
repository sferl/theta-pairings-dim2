# the miller function μ_{nP,nP}(Q)
def miller_double(nP, Q):
    l = nP._line_(nP, Q)
    nnP = nP+nP
    v = nnP._line_(-nnP, Q)
    return (nnP, l/v)

# the miller function μ_{nP,P}(Q)
def miller_add(nP, P, Q):
    l = nP._line_(P, Q)
    nPP = nP+P
    v = nPP._line_(-nPP, Q)
    return (nPP, l/v)

def miller_add(nP, P, Q):
    pass

def ladder0(n, P, Q):
    nP=P; r=1
    for bit in bin(n)[3:]:
        nP, s = miller_double(nP, Q)
        r=r*r*s
        if bit == "1":
            nP, s = miller_add(nP, P, Q)
            r=r*s
    return r


def tate0(n, P, Q, k=1, exp=ladder0):
    p=P.curve().base_ring().characteristic()
    assert (p**k-1)%n == 0
    d=(p**k-1)//n
    e = exp(n, P, Q)
    return e**d

if __name__ == "__main__" and "__file__" in globals():
    import time
    from sage.all import ZZ, GF, EllipticCurve, proof

    def testBN():
        p=16030569034403128277756688287498649515636838101184337499778392980116222246913;
        F=GF(p)
        E=EllipticCurve(F, [0,5])
        # E.j_invariant()=0, p=1 mod 4 so p splits in Q(i) and E ordinary
        r=E.cardinality() # this is prime: r.is_prime(): True
        #k=emb_degree(r,p) # 12
        k=12
        F2=GF(p**12, 't')
        E2=E.base_extend(F2)
        r2=E2.cardinality()

        #AC=(F2(6123150921968337422568796176365652246109592885288416056556751737102623758294), F2(1))
        #A24 = ((AC[0]/AC[1]+2)/4)
        #E2M=EllipticCurve(F2, [0,AC[0],0,1,0]) #this is a twist of E over Fp

        P=E.random_point()
        Q=(r2//r**2)*E2.random_point() # multiply by cofactor to get point of r-torsion
        print(f"- Test pairings on G1xG3, n={r}\n")
        #test_odd_degree_pairings(E2, r, r2//r**2, A24, k=k, points=[phi(E2(P)),phi(Q)])

        time0 = time.time()
        t1=tate0(r, E2(P), Q, k=k)
        print(f"Double And Add Miller ladder: {time.time() - time0:.5f}")
        time0 = time.time()

    print("########### Testing a BN curve ############")
    testBN()
