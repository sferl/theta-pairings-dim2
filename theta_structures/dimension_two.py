# Sage Imports
from sage.all import (
    cached_method,
    Integer,
    HyperellipticCurve,
    PolynomialRing,
    matrix
)
from sage.structure.element import get_coercion_model, RingElement
from utilities.batched_inversion import batched_inversion
from biextensions.morphism import Isogeny, TrivialChangeModel, Translation, LinearChangeModel


cm = get_coercion_model()


# ============================================ #
#     Class for Theta Structure (level-2?)     #
# ============================================ #


class ThetaStructure:
    """
    Class for the ThetaStructure, defined by its theta null point. This type
    represents the generic domain/codomain of the (2,2)-isogeny in the theta
    model.
    """

    def __init__(self, null_point):
        if not len(null_point) == 4:
            raise ValueError

        self._base_ring = cm.common_parent(*(c.parent() for c in null_point))
        self._point = ThetaPoint
        self._precomputation = None

        self._null_point = self._point(self, null_point)

    def null_point(self):
        """
        Return the null point of the given theta structure
        """
        return self._null_point

    def base_ring(self):
        """
        Return the base ring of the common parent of the coordinates of the null point
        """
        return self._base_ring

    def zero(self):
        """
        The additive identity is the theta null point
        """
        return self.null_point()

    def __repr__(self):
        return f"Theta structure over {self.base_ring()} with null point: {self.null_point()}"

    def coords(self):
        """
        Return the coordinates of the theta null point of the theta structure
        """
        return self.null_point().coords()

    def hadamard(self):
        """
        Compute the Hadamard transformation of the theta null point of the theta structure
        """
        return self.null_point().hadamard()

    def squared_theta(self):
        """
        Square the coefficients and then compute the Hadamard transformation of
        the theta null point of the theta structure
        """
        return self.null_point().squared_theta()

    def _arithmetic_precomputation(self):
        """
        Precompute 6 field elements used in arithmetic and isogeny computations
        """
        if self._precomputation is None:
            a, b, c, d = self.null_point().coords()

            # Technically this computes 4A^2, 4B^2, ...
            # but as we take quotients this doesnt matter
            # Cost: 4S
            AA, BB, CC, DD = self.squared_theta()

            # Precomputed constants for addition and doubling
            b_inv, c_inv, d_inv, BB_inv, CC_inv, DD_inv = batched_inversion(
                b, c, d, BB, CC, DD
            )

            y0 = a * b_inv
            z0 = a * c_inv
            t0 = a * d_inv

            Y0 = AA * BB_inv
            Z0 = AA * CC_inv
            T0 = AA * DD_inv

            self._precomputation = (y0, z0, t0, Y0, Z0, T0)
        return self._precomputation
    
    def translate_by(self, T):
        def translate_coords_by_index(P_coords, indices):
            a,b,c,d = P_coords
            i0,i1,j0,j1 = indices
            for i in range(i0):
                a,b,c,d = b,a,d,c
            for i in range(i1):
                a,b,c,d = c,d,a,b
            for j in range(j0):
                a,b,c,d = a,b,-c,-d
            for j in range(j1):
                a,b,c,d = a,-b,c,-d
            return a,b,c,d
        
        null_point = self.coords()
        for i in range(16):
            ind = (i//8, (i%8)//4, (i%4)//2, i%2)
            T_ind = self(translate_coords_by_index(null_point, ind))
            if T == T_ind:
                N = matrix(
                    [list(translate_coords_by_index((1,0,0,0), ind)),
                     list(translate_coords_by_index((0,1,0,0), ind)),
                     list(translate_coords_by_index((0,0,1,0), ind)),
                     list(translate_coords_by_index((0,0,0,1), ind))]
                )
                return Translation(domain=self, N=N)
        raise ValueError("The argument of this method should be a valid 2-torsion point.")

    @cached_method
    def rosenhain_from_theta(self):
        """
        From a theta null point structure, return the Rosenhain invariants

        NOTE: used for debugging
        """
        # Extract out the hadamard transform from the point class
        to_hadamard = self.zero().to_hadamard

        a, b, c, d = self.coords()
        A, C, B, D = to_hadamard(a, d, b, c)  # beware the weird order
        if A * B * C * D == 0:
            a, b, c, d = to_hadamard(a, b, c, d)
            A, C, B, D = to_hadamard(a, d, b, c)

        try:
            al = a * a + b * b + c * c + d * d
            be = 2 * (a * d + b * c)
            ga = 2 * (a * b + c * d)
            de = 2 * (a * c + b * d)
            CD_AB = C * D / (A * B)
            ep_phi = (1 + CD_AB) / (1 - CD_AB)
            lam = al * ga / (de * be)
            mu = ga / de * ep_phi
            nu = al / be * ep_phi

        except Exception as e:
            raise ValueError(f"Conversion to rosenhain failed because of: {e}")

        return lam, mu, nu

    @cached_method
    def hyperelliptic_from_theta(self):
        """
        Convert a theta null point structure to an hyperelliptic curve
        """

        lam, mu, nu = self.rosenhain_from_theta()
        if lam is None:
            raise ValueError("Could not compute Rosenhain roots from the null point")

        R = PolynomialRing(self.base_ring(), name="x")
        x = R.gens()[0]

        f_poly = x * (x - 1) * (x - lam) * (x - mu) * (x - nu)

        try:
            H = HyperellipticCurve(f_poly)
        except:
            raise ValueError("Converted curve is not a hyperelliptic curve")
        return H

    def __call__(self, P):
         # if coords == ThetaPoint(self, (0, 0, 0, 0)):
        if type(P) == tuple:
            coords = P
        elif isinstance(P.parent(), (ThetaStructure, AffineThetaStructure)):
            coords = P.coords()
        if coords == (0,0,0,0):
            raise ValueError("Cannot create a theta point with all zero coordinates")
        return self._point(self, coords)
    
    def to_affine(self):
        return TrivialChangeModel(self, AffineThetaStructure(self.coords()))


# ======================================== #
#     Class for Theta Point (level-2?)     #
# ======================================== #


class ThetaPoint:
    """
    A Theta Point in the level-2 Theta Structure is defined with four projective
    coordinates

    We cannot perform arbitrary arithmetic, but we can compute doubles and
    differential addition, which like x-only points on the Kummer line, allows
    for scalar multiplication
    """

    def __init__(self, parent, coords):
        if not isinstance(parent, ThetaStructure):
            raise ValueError

        self._parent = parent
        self._coords = tuple(coords)

        self._hadamard = None
        self._squared_theta = None

    def parent(self):
        """
        Return the parent of the element, of type ThetaStructure
        """
        return self._parent

    def theta(self):
        """
        Return the parent theta structure of this ThetaPoint"""
        return self.parent()

    def coords(self):
        """
        Return the projective coordinates of the ThetaPoint
        """
        return self._coords

    def is_zero(self):
        """
        An element is zero if it is equivalent to the null point of the parent
        ThetaStrcuture
        """
        return self == self.parent().zero()

    @staticmethod
    def to_hadamard(x_00, x_10, x_01, x_11):
        """
        Compute the Hadamard transformation of four coordinates, using recursive
        formula.
        """
        x_00, x_10 = (x_00 + x_10, x_00 - x_10)
        x_01, x_11 = (x_01 + x_11, x_01 - x_11)
        return x_00 + x_01, x_10 + x_11, x_00 - x_01, x_10 - x_11

    def hadamard(self):
        """
        Compute the Hadamard transformation of this element
        """
        if self._hadamard is None:
            self._hadamard = self.to_hadamard(*self.coords())
        return self._hadamard

    @staticmethod
    def to_squared_theta(x, y, z, t):
        """
        Square the coordinates and then compute the Hadamard transform of the
        input
        """
        return ThetaPoint.to_hadamard(x * x, y * y, z * z, t * t)

    def squared_theta(self):
        """
        Compute the Squared Theta transformation of this element
        which is the square operator followed by Hadamard.
        """
        if self._squared_theta is None:
            self._squared_theta = self.to_squared_theta(*self.coords())
        return self._squared_theta
    
    @staticmethod
    def _componentwise_multiply(P_coords, Q_coords):
        """
        Componentwise multiply two coordinate tuples
        """
        return ([a * b for a, b in zip(P_coords, Q_coords)])

    def double(self):
        """
        Computes [2]*self

        NOTE: Assumes that no coordinate is zero

        Cost: 8S 6M
        """
        # If a,b,c,d = 0, then the codomain of A->A/K_2 is a product of
        # elliptic curves with a non product theta structure.
        # Unless we are very unlucky, A/K_1 will not be in this case, so we
        # just need to Hadamard, double, and Hadamard inverse
        # If A,B,C,D=0 then the domain itself is a product of elliptic
        # curves with a non product theta structure. The Hadamard transform
        # will not change this, we need a symplectic change of variable
        # that puts us back in a product theta structure
        y0, z0, t0, Y0, Z0, T0 = self.parent()._arithmetic_precomputation()

        # Temp coordinates
        # Cost 8S 3M
        xp, yp, zp, tp = self.squared_theta()
        xp = xp**2
        yp = Y0 * yp**2
        zp = Z0 * zp**2
        tp = T0 * tp**2

        # Final coordinates
        # Cost 3M
        X, Y, Z, T = self.to_hadamard(xp, yp, zp, tp)
        X = X
        Y = y0 * Y
        Z = z0 * Z
        T = t0 * T

        coords = (X, Y, Z, T)
        return self._parent(coords)

    def diff_add(P, Q, PQ):
        """
        Given the theta points of P, Q and P-Q computes the theta point of
        P + Q.

        NOTE: Assumes that no coordinate is zero

        Cost: 8S 17M
        """
        # Extract out the precomputations
        Y0, Z0, T0 = P.parent()._arithmetic_precomputation()[-3:]

        # Transform with the Hadamard matrix and multiply
        # Cost: 8S 7M
        p1, p2, p3, p4 = P.squared_theta()
        q1, q2, q3, q4 = Q.squared_theta()

        xp = p1 * q1
        yp = Y0 * p2 * q2
        zp = Z0 * p3 * q3
        tp = T0 * p4 * q4

        # Final coordinates
        PQx, PQy, PQz, PQt = PQ.coords()

        # Note:
        # We replace the four divisions by
        # PQx, PQy, PQz, PQt by 10 multiplications
        # Cost: 10M
        PQxy = PQx * PQy
        PQzt = PQz * PQt

        X, Y, Z, T = P.to_hadamard(xp, yp, zp, tp)
        X = X * PQzt * PQy
        Y = Y * PQzt * PQx
        Z = Z * PQxy * PQt
        T = T * PQxy * PQz

        coords = (X, Y, Z, T)
        return P.parent()(coords)
    
    # we have P=self, Q, R, R+Q, P+R, P+Q and we compute P+Q+R
    def three_way_add(self, Q, R, QR, PR, PQ):
        # we define it only in the affine case
        raise NotImplementedError("This method is not implemented for projective theta structures. Use affine coordinates")
    
    # we have P=self, P+T, Q, Q+T and we compute P+Q, P+Q+T
    def compatible_add(self, PT, Q, QT):
        raise NotImplementedError("This method is not implemented yet.")

    def scale(self, n):
        """
        Scale all coordinates of the ThetaPoint by `n`
        """
        x, y, z, t = self.coords()
        if not isinstance(n, RingElement):
            raise ValueError(f"Cannot scale by element {n} of type {type(n)}")
        scaled_coords = (n * x, n * y, n * z, n * t)
        return self._parent(scaled_coords)
    
    def ratio(self, other):
        """
        Return the ratio of two theta points, lambda s.t. lambda P = Q
        """
        for i in range(4):
            if self.coords()[i] != 0:
                l = other.coords()[i]/self.coords()[i]
                assert self.scale(l).coords() == other.coords()
                return l

    def double_iter(self, m):
        """
        Compute [2^n] Self

        NOTE: Assumes that no coordinate is zero at any point during the doubling
        """
        if not isinstance(m, Integer):
            try:
                m = Integer(m)
            except:
                raise TypeError(f"Cannot coerce input scalar {m = } to an integer")

        if m.is_zero():
            return self.parent().zero()

        P1 = self
        for _ in range(m):
            P1 = P1.double()
        return P1

    def __mul__(self, m):
        """
        Uses Montgomery ladder to compute [m] Self

        NOTE: Assumes that no coordinate is zero at any point during the doubling
        """
        # Make sure we're multiplying by something value
        if not isinstance(m, (int, Integer)):
            try:
                m = Integer(m)
            except:
                raise TypeError(f"Cannot coerce input scalar {m = } to an integer")

        # If m is zero, return the null point
        if not m:
            return self.parent().zero()

        # We are with ±1 identified, so we take the absolute value of m
        m = abs(m)

        P0, P1 = self, self
        P2 = P1.double()
        # If we are multiplying by two, the chain stops here
        if m == 2:
            return P2

        # Montgomery double and add.
        for bit in bin(m)[3:]:
            Q = P2.diff_add(P1, P0)
            if bit == "1":
                P2 = P2.double()
                P1 = Q
            else:
                P1 = P1.double()
                P2 = Q

        return P1

    def __rmul__(self, m):
        return self * m

    def __imul__(self, m):
        self = self * m
        return self
    
    def translate_by(self, T):
        f = self.parent().translate_by(T)
        return f(self)
    
        # we assume we have P-Q; we could instead assume we have P+Q with some minor changes
    def full_ladder3(self, n, Q, PmQ):
        """
        From [self=P], [Q], [P-Q], return [nQ], [P+n Q]
        """
        P = self
        nQ=Q._parent.zero()
        nQQ=Q
        nQP=P

        if n == 0:
            return nQ, nQP

        if n < 0:
            PpQ = P.diff_add(Q, PmQ)
            return self.full_ladder3(abs(n), Q, PpQ)

        # Montgomery-ladder
        for bit in bin(n)[2:]:
            # on the first bit we add 0, so we could treat this case separately
            R = nQQ.diff_add(nQ, Q)
            if bit == "0":
                nQP=nQP.diff_add(nQ, P)
                nQ=nQ.double()
                nQQ=R
            else:
                nQP=nQP.diff_add(nQQ, PmQ)
                nQQ=nQQ.double()
                nQ=R
        return (nQ, nQP)

    # we assume here we have P+Q
    def full_ladder3_bis(self, n, Q, PQ):
        """
        From [self=P], [Q], [P+Q], return [nQ], [P+n Q]
        """
        if n == 0:
            return Q._parent.zero(), self
        if n == 1:
            return Q, PQ

        P = self
        nQ=Q
        nQQ=Q.double()
        nQQP=PQ.diff_add(Q, P)

        # if n == 2:
        #     return nQQ, nQQP

        if n < 0:
            PmQ = P.diff_add(Q, PQ)
            return self.full_ladder3_bis(abs(n), Q, PmQ)

        for bit in bin(n-1)[3:]:
            R = nQQ.diff_add(nQ, Q)
            if bit == "0":
                nQQP=nQQP.diff_add(nQ, PQ)
                nQ=nQ.double()
                nQQ=R
            else:
                nQQP=nQQP.diff_add(nQQ, P)
                nQQ=nQQ.double()
                nQ=R
        return (nQQ, nQQP)

    def ladder3(self, n, Q, PmQ):
        nQ, nQP=self.full_ladder3(n, Q, PmQ)
        return nQP

    def __eq__(self, other):
        """
        Check the quality of two ThetaPoints. Note that as this is a
        projective equality, we must be careful for when certain coefficients may
        be zero.
        """
        if not isinstance(other, ThetaPoint):
            return False

        a1, b1, c1, d1 = self.coords()
        a2, b2, c2, d2 = other.coords()

        if d1 != 0 or d2 != 0:
            return all([a1 * d2 == a2 * d1, b1 * d2 == b2 * d1, c1 * d2 == c2 * d1])
        elif c1 != 0 or c2 != 0:
            return all([a1 * c2 == a2 * c1, b1 * c2 == b2 * c1])
        elif b1 != 0 or b2 != 0:
            return a1 * b2 == a2 * b1
        else:
            return True

    def __repr__(self):
        return f"Theta point with coordinates: {self.coords()}"

class AffineThetaStructure(ThetaStructure):
    def __init__(self, constants, **kwds):
        super().__init__(constants, **kwds)
        self._point = AffineThetaPoint
        self._type = "AffineTheta"

    @cached_method
    def _arithmetic_precomputation(self):
        if self._precomputation is None:
            a, b, c, d = self.null_point().coords()

            # Technically this computes 4A^2, 4B^2, ...
            # but as we take quotients this doesnt matter
            # Cost: 4S
            AA, BB, CC, DD = self.squared_theta()

            # Precomputed constants for addition and doubling
            self._precomputation = batched_inversion(
                a, b, c, d, 4*AA, 4*BB, 4*CC, 4*DD
            )

        return self._precomputation

    def to_projective(self):
        return TrivialChangeModel(self, ThetaStructure(self.coords()))
    
    def translate_by(self, T):
        def translate_coords_by_index(P_coords, indices):
            a,b,c,d = P_coords
            i0,i1,j0,j1 = indices
            for i in range(i0):
                a,b,c,d = b,a,d,c
            for i in range(i1):
                a,b,c,d = c,d,a,b
            for j in range(j0):
                a,b,c,d = a,b,-c,-d
            for j in range(j1):
                a,b,c,d = a,-b,c,-d
            return a,b,c,d
        
        null_point = self.coords()
        proj = ThetaStructure(null_point)
        for i in range(16):
            ind = (i//8, (i%8)//4, (i%4)//2, i%2)
            T_ind = self(translate_coords_by_index(null_point, ind))
            if proj(T) == proj(T_ind):
                N = matrix(
                    [list(translate_coords_by_index((1,0,0,0), ind)),
                     list(translate_coords_by_index((0,1,0,0), ind)),
                     list(translate_coords_by_index((0,0,1,0), ind)),
                     list(translate_coords_by_index((0,0,0,1), ind))]
                )
                return Translation(domain=self, N=N)
        raise ValueError("The argument of this method should be a valid 2-torsion point.")
    
class AffineThetaPoint(ThetaPoint):
    def __eq__(self, other):
        """
        Check the quality of two ThetaPoints.
        """
        if isinstance(other, ThetaPoint):
            other = self._parent(other.coords())
        if not isinstance(other, AffineThetaPoint):
            return False
        return self.coords() == other.coords()
    
    def is_proj_eq(self, other):
        """
        Check the quality of two ThetaPoints. Note that as this is a
        projective equality, we must be careful for when certain coefficients may
        be zero.
        """
        if isinstance(other, ThetaPoint):
            other = self._parent(other.coords())
        if not isinstance(other, AffineThetaPoint):
            return False

        a1, b1, c1, d1 = self.coords()
        a2, b2, c2, d2 = other.coords()

        if d1 != 0 or d2 != 0:
            return all([a1 * d2 == a2 * d1, b1 * d2 == b2 * d1, c1 * d2 == c2 * d1])
        elif c1 != 0 or c2 != 0:
            return all([a1 * c2 == a2 * c1, b1 * c2 == b2 * c1])
        elif b1 != 0 or b2 != 0:
            return a1 * b2 == a2 * b1
        else:
            return True

    def diff_add(P, Q, PQ):
        """
        Given the theta points of P, Q and P-Q computes the theta point of
        P + Q.

        NOTE: Assumes that no coordinate is zero

        Cost: 8S 17M
        """
        # Extract out the precomputations
        Ainv, Binv, Cinv, Dinv = P.parent()._arithmetic_precomputation()[-4:]

        # Transform with the Hadamard matrix and multiply
        # Cost: 8S 7M
        p1, p2, p3, p4 = P.squared_theta()
        q1, q2, q3, q4 = Q.squared_theta()

        xp = Ainv * p1 * q1
        yp = Binv * p2 * q2
        zp = Cinv * p3 * q3
        tp = Dinv * p4 * q4

        # Final coordinates
        PQx_inv, PQy_inv, PQz_inv, PQt_inv = batched_inversion(*(PQ.coords()))

        X, Y, Z, T = P.to_hadamard(xp, yp, zp, tp)
        X = X * PQx_inv
        Y = Y * PQy_inv
        Z = Z * PQz_inv
        T = T * PQt_inv

        coords = (X, Y, Z, T)
        return P.parent()(coords)

    def double(self):
        """
        Computes [2]*self

        NOTE: Assumes that no coordinate is zero

        Cost: 8S 6M
        """
        # If a,b,c,d = 0, then the codomain of A->A/K_2 is a product of
        # elliptic curves with a non product theta structure.
        # Unless we are very unlucky, A/K_1 will not be in this case, so we
        # just need to Hadamard, double, and Hadamard inverse
        # If A,B,C,D=0 then the domain itself is a product of elliptic
        # curves with a non product theta structure. The Hadamard transform
        # will not change this, we need a symplectic change of variable
        # that puts us back in a product theta structure
        ainv, binv, cinv, dinv, Ainv, Binv, Cinv, Dinv = self.parent()._arithmetic_precomputation()

        # Temp coordinates
        # Cost 8S 3M
        xp, yp, zp, tp = self.squared_theta()
        xp = Ainv * xp**2
        yp = Binv * yp**2
        zp = Cinv * zp**2
        tp = Dinv * tp**2

        # Final coordinates
        # Cost 3M
        X, Y, Z, T = self.to_hadamard(xp, yp, zp, tp)
        X = ainv * X
        Y = binv * Y
        Z = cinv * Z
        T = dinv * T

        coords = (X, Y, Z, T)
        return self._parent(coords)
    
    def three_way_add(self, Q, R, QR, PR, PQ):
        S_Q_R = self._componentwise_multiply(Q.coords(), R.coords())
        S_Q_R = self.to_hadamard(*S_Q_R)

        S_0_QR = self._componentwise_multiply(self.parent().coords(), QR.coords())
        S_0_QR = self.to_hadamard(*S_0_QR)

        S_PQ_PR = self._componentwise_multiply(PQ.coords(), PR.coords())
        S_PQ_PR = self.to_hadamard(*S_PQ_PR)

        xQR_inv, yQR_inv, zQR_inv, tQR_inv, xP_inv, yP_inv, zP_inv, tP_inv = batched_inversion(*S_Q_R, *self.scale(4).coords())

        S_sums = self._componentwise_multiply(S_0_QR, S_PQ_PR)
        S_sums = self._componentwise_multiply(S_sums,
                                              (xQR_inv, yQR_inv, zQR_inv, tQR_inv))
        S_sums = self.to_hadamard(*S_sums)

        S_PQR = self._componentwise_multiply(S_sums, (xP_inv, yP_inv, zP_inv, tP_inv))

        return self._parent(S_PQR)
    
    # we have P=self, P+T, Q, Q+T and we compute P+Q, P+Q+T
    def compatible_add(self, PT, Q, QT):
        raise NotImplementedError("This method is not implemented yet.")
        
    def translate_by(self, T):
        f = self.parent().translate_by(T)
        return f(self)