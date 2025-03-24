from biextensions.biextension import Biextension
from sage.all import proof, next_prime, supersingular_j
from sage.schemes.elliptic_curves.constructor import EllipticCurve, EllipticCurve_from_j
from sage.rings.finite_rings.finite_field_constructor import GF
from sage.rings.integer_ring import ZZ
from sage.groups.generic import order_from_multiple
from utilities.discrete_log import discrete_log_pari
from theta_structures.couple_point import CouplePoint
from theta_isogenies.gluing_isogeny import GluingThetaIsogeny
from theta_isogenies.isogeny import ThetaIsogeny

proof.all(False)

##################################
# Utilities for testing
##################################

def ss_torsion_basis(E, n, cofactor=1):
    """
    Compute a basis for the n-torsion of E.
    where E is a supersingular elliptic curve
    with rational n*cofactor torsion.
    """

    def point_order_n(E, n, cofactor):
        """
        Compute a random point P of order n
        """
        while True:
            P = cofactor * E.random_point()
            P.set_order(multiple=n)
            if P.order() == n:
                return P

    # We can pick any point P of order n and then find some
    # Q such that e(P, Q) has full order
    P = point_order_n(E, n, cofactor)
    while True:
        Q = point_order_n(E, n, cofactor)
        ePQ = P.weil_pairing(Q, n)
        o = order_from_multiple(ePQ, n, operation="*")
        if o == n:
            return P, Q
        
def gen_params(e):
    """
    Generate a prime p of the form (2**e * r * cofactor - 1) where r is a prime of e bits
    """
    cofactor = ZZ(1)
    r = next_prime(ZZ(1 << (e-1)))
    while True:
        p = ZZ(cofactor << e) * r - 1
        if p.is_prime():
            Fp2 = GF(p**2, name='i', modulus=[1,0,1])
            return p, r, cofactor, Fp2
        cofactor += 1


# Q is the n-torsion point
def compute_tate_pairings(n, P, Q, PQ, k=1, d=None, exp_function=None, zero=None, scale = False):
    # n is the order of the point Q
    # P is any point on a Kummer variety K
    # Q is an n-torsion point on K
    # PQ is P+Q
    # k is the embedding degree wrt the base field and n
    # d is the final exponentiation exponent
    # exp_function is g -> g^n in the second component
    # zero is the null point
    # scale is a boolean telling whether to apply random scaling to the points
    # we keep returning t1 for compatibility with the original repository. TODO: fix this
    
    # NOTE: with respect to compute_tate_pairings on Damien's repo, this function takes in an additional parameter PQ and doesn't take in the parent variety K

    # NOTE: for even pairing, we need to start with an element g_P,Q which is a
    # tensorial square of h_P,Q in level 1; hence if we scale we should
    # scale by an element which is globally a square over the base ring
    
    R=P.parent().base_ring()
    p = R.characteristic()
    if scale:
        # random scalings to test the algorithm
        lambdaP=R.random_element()
        lambdaQ=R.random_element()
        lambdaPQ=R.random_element()
        P=P.scale(lambdaP)
        Q=Q.scale(lambdaQ)
        PQ=PQ.scale(lambdaPQ)
        print([t.is_square() for t in [lambdaP, lambdaQ, lambdaPQ]])

    g=Biextension(P, Q, PQ, zero=zero)

    # t1=Q.tate_pairing(P, n, k=k)
    # TODO:
    # NEED A REFERENCE Tate PAIRING FOR HYPERELLIPTIC CURVES
    t1 = 1


    (t2, t2a)=g.tate_pairing(n, k=k, d=d, exp_function=exp_function)
    # the reduced tate pairing is also given by using the exponentiation by p^k-1 in the biextension
    t2bis = g.non_reduced_tate_pairing(p**k-1, exp_function=exp_function)
    return t1, t2, t2a, t2bis

def compute_weil_pairings(n, P, Q, PQ, k=1, exp_function=None, zero=None, scale=False):
    # n is the order of the point Q
    # P is any point on a Kummer variety K
    # Q is an n-torsion point on K
    # PQ = P + Q
    # k is the embedding degree wrt the base field and n
    # d is the final exponentiation exponent
    # exp_function is g -> g^n in the second component
    # zero is the null point
    # scale is a boolean telling whether to apply random scaling to the points
    # we keep returning t1 for compatibility with the original repository. TODO: fix this

    if scale:
        # random scalings to test the algorithm
        R=P.parent().base_ring()

        lambdaP=R.random_element()
        lambdaQ=R.random_element()
        lambdaPQ=R.random_element()

        P=P.scale(lambdaP)
        Q=Q.scale(lambdaQ)
        PQ=PQ.scale(lambdaPQ)

    # t1=Q.weil_pairing(P, n)
    # TODO: NEED A REFERENCE WEIL PAIRING FOR HYPERELLIPTIC CURVES
    t1 = 1

    t2=Biextension(P, Q, PQ, zero=zero).weil_pairing(n, exp_function=exp_function)
    return t1, t2

def isogeny_chain_2pow_from_product(P1xP2, Q1xQ2, nsteps):
    # P1xP2 = CouplePoint(P1, P2)
    # Q1xQ2 = CouplePoint(Q1, Q2)
    assert order_from_multiple(P1xP2, ZZ(1 << (nsteps + 2)), operation="*") == ZZ(1 << (nsteps + 2))
    assert order_from_multiple(Q1xQ2, ZZ(1 << (nsteps + 2)), operation="*") == ZZ(1 << (nsteps + 2))

    E1, E2 = P1xP2.curves()
    assert (E1, E2) == Q1xQ2.curves()
    
    TP = P1xP2
    TQ = Q1xQ2
    for j in nsteps:
        kerP = ZZ(1 << (nsteps + 1 - j)) * TP
        kerQ = ZZ(1 << (nsteps + 1 - j)) * TQ
        assert order_from_multiple(kerP, ZZ(8), operation="*") == ZZ(8)
        if j == 1:
            # gluing step
            phi = GluingThetaIsogeny(kerP, kerQ)
        else:
            phi = ThetaIsogeny(Kum, kerP, kerQ) * phi
        Kum = phi.codomain()
        TP, TQ = phi(TP), phi(TQ)
    return phi, Kum

def generate_kummer(F, e):
    # Generate a Kummer variety K (both as a projective ThetaStructure Kum_proj and as an AffineThetaStructure Kum),
    # a (2,2) gluing isogeny phi from a product E x E' of supersingular elliptic curves to K,
    # a prime p of the form p = 2^e * r * f - 1, with e >= 2, r a large prime, f an integer cofactor
    # two points P, Q of order r on E x E' (so that phi(P), phi(Q) are points of order r on K), and a random point R on E x E'
    # returns: Kum_proj, Kum, phi, P, Q, R, p, e, r, f

    # take two random-ish supersingular elliptic curves
    p = F.characteristic()
    E0 = EllipticCurve_from_j(supersingular_j(F))
    P0, Q0 = ss_torsion_basis(E0, ZZ(1 << e), cofactor=((p+1) >> e))
    E1 = E0.isogeny(P0, algorithm="factored").codomain()
    E2 = E0.isogeny(Q0, algorithm="factored").codomain()
    
    # now walk nsteps steps, with nsteps in [1, e-2], in the (2,2)-isogeny graph to get a random-ish Kummer surface
    P1, Q1 = ss_torsion_basis(E1, ZZ(1 << e), cofactor=((p+1) >> e))
    P2, Q2 = ss_torsion_basis(E2, ZZ(1 << e), cofactor=((p+1) >> e))

    # find isotropic kernel subgroup <[2**(e-2-nsteps)](P1,P2), [2**(e-2-nsteps)](Q1,Q2)> :
    # i.e., we do nsteps in [1, e-2], and adjust the kernel to have order 2**(nsteps + 2)
    # in principle, we could do in total e steps of the isogeny to get a more random Kummer variety; we only do one for now.
    isotropic_factor = -discrete_log_pari(P2.weil_pairing(Q2, 2**e), P1.weil_pairing(Q1, 2**e), 2**e)
    P1xP2 = CouplePoint(isotropic_factor*P1, P2)
    Q1xQ2 = CouplePoint(Q1, Q2)
    nsteps = e-2
    assert 1 <= nsteps <= e-2
    TP, TQ = P1xP2 * ZZ(1 << (e - 2 - nsteps)), Q1xQ2 * ZZ(1 << (e - 2 - nsteps))
    phi, Kum_proj = isogeny_chain_2pow_from_product(TP, TQ, nsteps)

    Kum_aff = Kum_proj.to_affine().codomain()

    # TODO: get a basis of the r-torsion of the Kummer variety
    ## TODO: write random_element in theta coordinates

    # # get a basis of the r-torsion of the Kummer variety
    # P1, Q1 = torsion_basis(E1, r, cofactor=(p+1)//r)
    # P2, Q2 = torsion_basis(E2, r, cofactor=(p+1)//r)
    
    # P1_rand = E1.random_point()
    # P2_rand = E2.random_point()

    # P = CouplePoint(P1, P2)
    # Q = CouplePoint(Q1, Q2)
    # # check P, Q are indeed of order r
    # assert(Q * r == CouplePoint(E1.zero(), E2.zero()))
    # assert(P * r == CouplePoint(E1.zero(), E2.zero()))
    # R = CouplePoint(P1_rand, P2_rand)
    
    return Kum_proj, Kum_aff, phi

def test_pairings():
    # Compute a product of supersingular curves, go to an isogenous Jacobian and take its Kummer
    e = 50
    p, r, f, Fp2 = gen_params(e)
    Kum_proj, Kum_aff, phi = generate_kummer(p, f, Fp2)
    E1, E2 = phi.domain.E1, phi.domain.E2

    # TODO: generate basis of the r torsion on Kum_proj and corresponding points(?) on Kum_aff
    # TODO: generate a random point on Kum_proj and Kum_aff for Tate pairing
    # gives an the affine Kummer point phi(R) for a point R on E x E'
    def to_aff(T):
        return Kum_aff(phi(T).coords())
    
    # PK, QK are points of order r on Kum
    # RK is a random point on Kum
    # RQK = RK + QK, PQK = PK + QK
    PK = to_aff(P)
    QK = to_aff(Q)
    RK = to_aff(R)
    RQK = to_aff(R+Q)
    PQK = to_aff(P+Q)
    assert Kum_proj(QK.coords()) != Kum_proj.zero()
    assert Kum_proj((QK * r).coords()) == Kum_proj.zero()

    # call functionalities from biextensions.biextension
    def exp_function(g, n):
        return g.ladder(n)
    def exp_function_bis(g, n):
        _, g2=g.full_ladder_bis(n)
        return g2
    def fast_ladder(g, n):
        return g.fast_ladder(n)

    # compute the pairings in theta coordinates and check they give the expected result
    print("- Test pairings in theta via ladder3_bis")
    
    t1 = R.tate_pairing(Q, r) ** (p**2//r)
    _, t2, t2bis = compute_tate_pairings(r, RK, QK, RQK, k=2, scale=True)
    assert t2 == t2bis
    assert t1**4 == t2 or t1**(-4) == t2
    # the sign ambiguity is due to the fact that we work with Kummer points
    # a factor 2 at the exponent is given by the algorithm itself: the biextension pairing with level 2 coordinates is the square of the usual pairing
    # another factor 2 is given by the fact that RK = phi(R), QK = phi(Q) are pushed through a 2-isogeny
    
    _, t4, _ = compute_tate_pairings(r, RK, to_aff(Q+Q), to_aff(R+Q+Q), k=2, scale=True)
    assert t4 == t2*t2
    t5, t6, _ = compute_tate_pairings(r, to_aff(R+R), QK, to_aff(R+R+Q), k=2, scale=True)
    assert t6 == t2*t2

    # Weil pairing in theta
    _, w2 = compute_weil_pairings(r, PK, QK, PQK, scale=True)
    w1 = P.weil_pairing(Q, r)
    assert w1**4 == w2 or w1**(-4) == w2
    

    print("- Test pairings in theta via ladder3")
    _, t2, t2bis = compute_tate_pairings(r, RK, QK, RQK, k=2, exp_function = fast_ladder, scale=True)
    assert t2 == t2bis
    assert t1**4 == t2 or t1**(-4) == t2

    print("- Test pairings in theta via the biextension ladder")
    _, t2, t2bis = compute_tate_pairings(r, RK, QK, RQK, k=2, exp_function=exp_function, scale=True)
    assert t2 == t2bis
    assert t1**4 == t2 or t1**(-4) == t2

    print("- Test pairings in theta via the biextension ladder bis")
    _, t2, t2bis = compute_tate_pairings(r, RK, QK, RQK, k=2, exp_function = exp_function_bis, scale=True)
    assert t2 == t2bis
    assert t1**4 == t2 or t1**(-4) == t2

def generate_kummer_even():
    # Generate a Kummer variety K (both as a projective ThetaStructure Kum_proj and as an AffineThetaStructure Kum),
    # a (2,2) gluing isogeny phi from a product E x E' of supersingular elliptic curves to K,
    # a prime p of the form p = 2^e * 3^2 * r * f - 1, with e >= 2, r a large prime, f an integer cofactor
    # two points P_4r, Q_4r of order 4r on E x E' lying above ker(phi) (so that phi(P_4r), phi(Q_4r) are points of order 2r on K), and a random point R on E x E'
    # returns: Kum_proj, Kum, phi, P_4r, Q_4r, R, p, e, r, f

    # generate prime p
    e = 10 # e >= 2, so that p = 3 mod 4
    r = next_prime(2**100)
    f = 1
    p = 2**e * 3**2 * r * f - 1
    while not p.is_prime():
        f = f + 1
        p = 2**e * 3**2 *  r * f - 1

    F = GF(p**2, name='i', modulus=[1,0,1])
    E0 = EllipticCurve(F, [1,0])
    
    # generate E x E' from E0, for a little randomization
    P0, Q0, = E0.torsion_basis(2**e)
    P3, Q3 = E0.torsion_basis(3)
    phi1 = E0.isogeny(P3)
    phi2 = E0.isogeny(Q3)
    E1 = phi1.codomain()
    E2 = phi2.codomain()
    P1, Q1 = phi1(P0), phi1(Q0)
    P2, Q2 = phi2(P0), phi2(Q0)

    # find isotropic kernel subgroup generated by K8_1, K8_2: it will be the kernel of a gluing isogeny ExE' -> K
    # in principle, we could do in total e steps of the isogeny to get a more random Kummer variety; we only do one for now.
    isotropic_factor = -discrete_log_pari(P2.weil_pairing(Q2, 2**e), P1.weil_pairing(Q1, 2**e), 2**e)
    K8_1 = (2**(e-3)) * CouplePoint(isotropic_factor*P1, P2)
    K8_2 = (2**(e-3)) * CouplePoint(Q1, Q2)
    phi = GluingThetaIsogeny(K8_1, K8_2)
    
    Kum_proj = phi.codomain()
    Kum = Kum_proj.to_affine().codomain()


    # get a basis of the 2r-torsion of the Kummer variety
    P1_r, Q1_r = gen_torsion_basis(E1, r)
    P2_r, Q2_r = gen_torsion_basis(E2, r)
    
    P1_rand = E1.random_point()
    P2_rand = E2.random_point()

    P_4r = CouplePoint(P1_r, P2_r) + K8_1 * 2
    Q_4r= CouplePoint(Q1_r, Q2_r) + K8_2 * 2
    R = CouplePoint(P1_rand, P2_rand)
    
    # check P, Q are indeed of order 4r
    assert(Q_4r * 4* r == CouplePoint(E1.zero(), E2.zero()))
    assert(P_4r * 4 * r == CouplePoint(E1.zero(), E2.zero()))

    # just check that phi is indeed an isogeny
    assert phi(P_4r * 2) == phi(P_4r) * 2
    assert phi(R * 19) == phi(R) * 19
    assert phi(Q_4r)*4*r == phi(Q_4r*4*r)
    assert phi(Q_4r)*4*r == Kum_proj.zero()

    return Kum_proj, Kum, phi, P_4r, Q_4r, R, p, e, r, f

def test_pairings_even():
    # Compute a product supersingular curves ExE', go to an isogenous Jacobian and take its Kummer
    Kum_proj, Kum, phi, P_4r, Q_4r, R, p, e, r, f = generate_kummer_even()
    P = P_4r * 2
    Q = Q_4r * 2

    # gives an the affine Kummer point phi(R) for a point R on E x E'
    def to_aff(T):
        return Kum(phi(T).coords())
    
    # P, Q are points of order 2r on ExE'
    assert P * 2 *r == P * 0
    assert Q * 2 *r == P * 0
    assert phi(P) * r == phi(P * 2* r)

    # PK_2r, QK_2r are points of order 2r on Kum
    # PK, QK are their doubles, of order r on Kum
    # RK is a random point on Kum
    # RQK = RK + QK, PQK = PK + QK   
    RK = to_aff(R)
    RQK = to_aff(R+Q)
    PQK = to_aff(P+Q)
    PK_2r = to_aff(P_4r)
    QK_2r = to_aff(Q_4r)
    PK = PK_2r * 2
    QK = QK_2r * 2
    
    assert Kum_proj((PK * r).coords()) == Kum_proj.zero()
    assert Kum_proj((QK * r).coords()) == Kum_proj.zero()

    # call functionalities from biextensions.biextension
    def exp_function(g, n):
        return g.ladder(n)
    def exp_function_bis(g, n):
        _, g2=g.full_ladder_bis(n)
        return g2
    def fast_ladder(g, n):
        return g.fast_ladder(n)

    # compute the pairings in theta coordinates and check they give the expected result
    print("- Test pairings in theta via ladder3_bis")
    t1 = R.tate_pairing(Q*2, r) ** (p**2//r)
    """ 
    _, t2, t2bis = compute_tate_pairings(r, RK, QK, RQK, k=2, scale=True)
    assert t2 == t2bis
    assert t2 == t1**2 or t2 == t1**(-2)
    # the sign ambiguity is due to the fact that we work with Kummer points
    # a factor 2 at the RHS exponent is given by the algorithm itself: the biextension pairing with level 2 coordinates is the square of the usual pairing
    # another factor 2 is given by the fact that RK = phi(R), QK = phi(Q) are pushed through a 2-isogeny
    # a factor 2 at the LHS exponent is given by the fact t1 is computed on Q*2, while t2 is computed on phi(Q).
    # in the end, we're left with exponent 2 at the RHS.

    _, t4, _ = compute_tate_pairings(r, RK, to_aff(Q+Q), to_aff(R+Q+Q), k=2, scale=True)
    assert t4 == t2*t2
    t5, t6, _ = compute_tate_pairings(r, to_aff(R+R), QK, to_aff(R+R+Q), k=2, scale=True)
    assert t6 == t2*t2

    # Weil pairing in theta
    _, w2 = compute_weil_pairings(r, PK, QK, PQK, scale=True)
    w1 = (P*2).weil_pairing(Q*2, r)
    assert w2 == w1**(-1) or w2 == w1
    # here another exponent 2 at the LHS is given by the scalar 2 multiplying P and cancels out with the exponent 2 at the RHS

    print("- Test pairings in theta via ladder3")
    _, t2, t2bis = compute_tate_pairings(r, RK, QK, RQK, k=2, exp_function = fast_ladder, scale=True)
    assert t2 == t2bis
    assert t2 == t1**2 or t2 == t1**(-2)

    print("- Test pairings in theta via the biextension ladder")
    _, t2, t2bis = compute_tate_pairings(r, RK, QK, RQK, k=2, exp_function=exp_function, scale=True)
    assert t2 == t2bis
    assert t2 == t1**2 or t2 == t1**(-2)

    print("- Test pairings in theta via the biextension ladder bis")
    _, t2, t2bis = compute_tate_pairings(r, RK, QK, RQK, k=2, exp_function = exp_function_bis, scale=True)
    assert t2 == t2bis
    assert t2 == t1**2 or t2 == t1**(-2)
    """
    # the following algorithm computes an even pairing: the biextension algorithm here is the usual pairing, not its square
    print("- Test even pairings in theta via ladder3_bis")
    
    t1 = R.tate_pairing(Q_4r*2, 2*r) ** (p**2//(2*r))
    _, t2, t2a, t2bis = compute_tate_pairings(r*2, RK, QK_2r, to_aff(R+Q_4r), k=2, scale=True)
    assert t2 == t1 or t2 == t1**(-1) or t2a == t1 or t2a == t1**(-1)
    # TODO: doesn't work when scaling with a nonsquare is applied
    """
    _, w2 = compute_weil_pairings(r*2, PK_2r, QK_2r, to_aff(P_4r + Q_4r), k=2)
    w1 = (P_4r*2).weil_pairing(Q_4r*2, 2*r)
    assert w2**(-2) == w1 or w2**2 == w1
    """
# test_pairings()
test_pairings_even()