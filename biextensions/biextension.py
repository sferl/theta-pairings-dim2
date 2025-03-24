from sage.all import Integer

class Biextension:
    @classmethod
    def diagonal(self, P, zero = None):
        return self(P, P, P.double(), zero=zero)

    # A biextension element is represented by
    # [0~, P~; Q~, P+Q~]
    # the biextension exponentiation implemented is the one which gives
    # [0~, P~, nQ~, P+nQ~]
    def __init__(self, P, Q, PQ, zero = None):
        assert Q.parent() == P.parent()
        assert PQ.parent() == P.parent()
        if not zero is None:
            # we need to check that 'zero' is a valid affine lift of the
            # neutral point
            assert zero == P.parent().zero()
        self._kummer = P.parent()
        self._P = P
        self._Q = Q
        self._PQ = PQ
        self._zero = zero

    def __repr__(self):
        return f"Biextension element: P={self.P().coords()}, Q={self.Q().coords()}, P+Q={self.PQ().coords()}"

    def parent(self):
        return self._kummer

    def zero(self):
        if self._zero is None:
            return self.parent().zero()
        return self._zero
    def P(self):
        return self._P
    def Q(self):
        return self._Q
    def PQ(self):
        return self._PQ

    def ratios(self, other):
        l0 = self.zero().ratio(other.zero())
        lP = self.P().ratio(other.P())
        lQ = self.Q().ratio(other.Q())
        lPQ = self.PQ().ratio(other.PQ())
        return (l0, lP, lQ, lPQ)

    def ratio(self, other):
        l0, lP, lQ, lPQ = self.ratios(other)
        return lPQ*l0 / (lP*lQ)

    def __eq__(self, other):
        if self.zero() != other.zero():
            return False
        if self.P() != other.P():
            return False
        if self.Q() != other.Q():
            return False
        if self.PQ() != other.PQ():
            return False
        l0, lP, lQ, lPQ = self.ratios(other)
        return l0*lPQ == lP*lQ

    def swap(self):
        return self.__class__(self.Q(), self.P(), self.PQ(), zero = self._zero)

    def new_element(self, Q, PQ):
        return self.__class__(self.P(), Q, PQ, zero = self._zero)

    #(0,0)
    def full_neutral(self):
        K=self.parent()
        zero = K.zero()
        return self.__class__(zero, zero, zero, zero = zero)

    #(P,0)
    def neutral(self):
        return self.new_element(self.zero(), self.P())

    # our arithmetic is done on the second argument
    # but we could have different affine lifts of 0, P on self and other,
    # so we need to adjust for that
    # in practice if we start with self, and do arithmetic on it, the
    # adjustment will be trivial
    # the returned adjustment 'a' is such that other.PQ().scale(a)
    # corresponds to the same P~
    def adjust(self, other):
        assert self.P() == other.P()
        l0 = self.zero().ratio(other.zero())
        lP = self.P().ratio(other.P())
        return lP/l0

    # by convention, we do our additions on the second argument Q
    def add(self, other):
        raise NotImplementedError("Biextension addition not implemented yet.")
        P=self.P()
        adj = self.adjust(other)
        Q1Q2, PQ1Q2_ = P.compatible_add(self.PQ(), other.Q(), other.PQ())
        # PQ1Q2_ is not correctly normalised
        PQ1Q2 = P.three_way_add(Q1, Q2, Q1Q2, PQ2, PQ1)
        assert PQ1Q2 == PQ1Q2_
        return self.new_element(Q1Q2, PQ1Q2.scale(adj))

    def double(self):
        QQ=self.Q().double()
        PQQ=self.PQ().diff_add(self.Q(), self.P())
        return self.new_element(QQ, PQQ)

    # the opposite is given by [P~, P+Q~; 0~ Q~]=[0~, Q~; -P~, P-Q~]
    def opposite(self):
        PmQ=self.P().diff_add(self.Q(), self.PQ())
        return self.new_element(self.Q(), PmQ)

    def __add__(self, other):
        if self == other:
            return self.double()
        return self.add(other)

    def __sub__(self, other):
        opp = other.opposite()
        return self.add(opp)

    # double and add algorithm
    # this computes nQ, P+nQ
    def mult(self, n):
        if n == 0:
            return self.neutral()
        if n < 0:
            opp = self.opposite()
            return opp.mult(abs(n))
        g = self
        for bit in bin(n)[2:]:
            g = self.double()
            if bit == "1":
                g = g.add(self)
        return g

    def __mul__(self, n):
        if not isinstance(n, (int, Integer)):
            try:
                n = Integer(n)
            except:
                raise TypeError(f"Cannot coerce input scalar {n = } to an integer")
        return self.mult(n)

    def __rmul__(self, n):
        return self * n

    # we are given [0~, P~; Q1~ P+Q1~] and [0~, P~; Q2~ P+Q2~] and
    # [0~, P~; Q1-Q2~, P+Q1-Q2~]
    # and we compute [0~, P~; Q1+Q2~ P+Q1+Q2~]
    def diff_add(self, other, diff):
        # we assume that we have the same 0~, P~ on self, other, diff
        # and that diff = self other^-1
        assert self.P() == (other.P())
        assert self.P() == (diff.P())

        Q12=self.Q().diff_add(other.Q(), diff.Q())
        PQ12=self.PQ().diff_add(other.Q(), diff.PQ())
        g1=self.new_element(Q12, PQ12)

        # sanity checks
        oppdiff = diff.opposite()
        Q12bis=other.Q().diff_add(self.Q(), oppdiff.Q())
        PQ12bis=other.PQ().diff_add(self.Q(), oppdiff.PQ())
        g2=self.new_element(Q12, PQ12)
        assert g1 == g2

        # note that we only used Q1, Q2, Q1mQ2, PQ1, PQ1mQ2
        # and that PQ2 is only there for the sanity check
        return g1

    # returns n.g, (n+1).g
    def full_ladder(self, n):
        if n < 0:
            opp = self.opposite()
            return opp.biext_full_ladder(abs(n))

        g = self
        g1 = self.neutral()
        g2 = self

        if n == 0:
            return g1, g2

        for bit in bin(n)[2:]:
            # on the first bit we add 0, so we could treat this case separatly
            if bit == "0":
                g2 = g2.diff_add(g1, g)
                g1 = g1.double()
            else:
                g1 = g2.diff_add(g1, g)
                g2 = g2.double()

        return g1, g2

    #like full_ladder, but start the ladder with g, 2g
    #and return (n-1).g, n.g rather than n.g, (n+1).g
    #(cf the comment in 'kummer')
    def full_ladder_bis(self, n):
        if n < 0:
            opp = self.opposite()
            return opp.biext_full_ladder(abs(n))

        g = self
        g1 = self
        g2 = self.double()

        if n == 0:
            return self, self.neutral()
        if n == 1:
            return self.neutral(), self
        #if n == 2:
        #    return g1, g2

        for bit in bin(n-1)[3:]:
            if bit == "0":
                g2 = g2.diff_add(g1, g)
                g1 = g1.double()
            else:
                g1 = g2.diff_add(g1, g)
                g2 = g2.double()

        return g1, g2

    def ladder(self, n):
        g1, g2 = self.full_ladder(n)
        return g1

    def translate_by(self, T):
        QT=self.Q().translate_by(T)
        PQT=self.PQ().translate_by(T)
        return self.new_element(QT, PQT)

    # in the ladder, we have g^n: nQ~, P+nQ~, g^{n+1}: (n+1)Q~, P+(n+1)Q~
    # but we can represent g^{n+1} only from (n+1)Q~;
    # P+(n+1)Q~ is (implicitly) given by the three way add of P, nQ, Q; P+nQ, P+Q, (n+1)Q
    # this 'compact' representation is exactly what 'full_ladder3' computes in the Kummer model; this saves one diff_add by operation
    def fast_ladder(self, n):
        # ladder3 assumes we have P, Q, P-Q
        # but the biext assumes we have P, Q, P+Q
        nQ, nQP = self.P().full_ladder3(-n, self.Q(), self.PQ())
        return self.new_element(nQ, nQP)

    def fast_ladder_bis(self, n):
        nQ, nQP = self.P().full_ladder3_bis(n, self.Q(), self.PQ())
        return self.new_element(nQ, nQP)

    # here self.Q() is a point of n-torsion
    def non_reduced_tate_pairing(self, n, exp_function=None):
        if exp_function is None:
            def exp_function(g, n):
                return g.fast_ladder_bis(n)
        g = exp_function(self, n)
        return self.neutral().ratio(g)

    # special case when n is even: we can use the action of G(20E) to
    # compute the Tate pairing on (0E) rather than 2(0E)
    def even_non_reduced_tate_pairing(self, n, exp_function=None):
        assert n%2==0
        m=n//2
        if exp_function is None:
            def exp_function(g, n):
                return g.fast_ladder_bis(n)
        g = exp_function(self, m)
        Q = g.Q()
        gT = g.translate_by(Q)
        return self.neutral().ratio(gT)

    def tate_pairing(self, n, k=1, d=None, exp_function=None):
        if d is None:
            p=self.P().parent().base_ring().characteristic()
            d=(p**k-1)/n
        if n%2==0:
            r=self.even_non_reduced_tate_pairing(n, exp_function=exp_function)
        else:
            r=self.non_reduced_tate_pairing(n, exp_function=exp_function)
        return r**d

    def weil_pairing(self, n, exp_function=None):
        if n%2==0:
            r1=self.even_non_reduced_tate_pairing(n, exp_function=exp_function)
            r2=self.swap().even_non_reduced_tate_pairing(n, exp_function=exp_function)
        else:
            r1=self.non_reduced_tate_pairing(n, exp_function=exp_function)
            r2=self.swap().non_reduced_tate_pairing(n, exp_function=exp_function)
        return r1/r2
