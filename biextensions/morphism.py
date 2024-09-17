from sage.all import Matrix

class Morphism:
    def __init__(self, domain = None, codomain = None, image=None, inverse_image = None):
        self._domain = domain
        self._codomain = codomain
        self._image = image
        self._inverse_image = inverse_image
        pass

    def domain(self):
        return self._domain

    def codomain(self):
        return self._codomain

    def image(self, P):
        return self._image(P)

    def inverse_image(self, P):
        return self._inverse_image(P)

    def dual(self):
        return self.__class__(domain = self.codomain(), codomain = self.domain(), image = self.inverse_image, inverse_image = self.image)

    def inverse(self):
        return self.dual()

    def __repr__(self):
        return f'Morphism {self.domain()} -> {self.codomain()}'

    def __call__(self, P):
        #print(f"Domain: {self.domain()}, codomain: {self.codomain()}, P: {P}")
        if not isinstance(P, self.domain().__class__):
            P=self.domain()(P)
        Q = self.image(P)
        if not isinstance(Q, self.codomain().__class__):
            Q=self.codomain()(Q)
        return Q

    def compose(self, f2):
        return ComposedMorphism(f1=self, f2=f2)

    def __mul__(self, psi):
        return self.compose(psi)

    def __rmul__(self, psi):
        return psi.compose(self)

# an isomorphism is supposed to have an inverse image
# note that a Morphism can have one too without being an isomorphism, think
# about the reduction from E to its Kummer line, the 'inverse' is lifting
# back to E
class Isomorphism(Morphism):
    @classmethod
    def isomorphism(cls, E1, E2):
        phi=E1.isomorphism_to(E2)
        phiinv=phi**(-1)
        def image(P):
            return phi(P)
        def inverse_image(P):
            return phiinv(P)
        return cls(domain=E1, codomain=E2, image=image, inverse_image=inverse_image)

class ComposedMorphism(Morphism):
    def __init__(self, f1 = None, f2 = None):
        self._domain = f1.domain()
        self._codomain = f2.codomain()
        self._f1 = f1
        self._f2 = f2

    def first(self):
        return self._f1
    def second(self):
        return self._f2

    def f1(self, P):
        return self.first()(P)
    def f2(self, P):
        return self.second()(P)

    def dual(self):
        f1dual=self._f1.dual()
        f2dual=self._f2.dual()
        return ComposedMorphism(f2dual, f1dual)

    def image(self, P):
        P1 = self.f1(P)
        P2 = self.f2(P1)
        return P2

    def inverse_image(self, P):
        P2 = self.second().inverse_image(P)
        P1 = self.first().inverse_image(P2)
        return P1

class LinearIsomorphism(Isomorphism):
    def __init__(self, domain = None, codomain = None, N = None, **kwds):
        super().__init__(domain=domain, codomain=codomain,**kwds)
        self._N = N

    def N(self):
        return self._N

    def dual(self):
        if self.N() is None:
            raise ValueError("Cannot compute the dual without the matrix N")
        N_inverse = self.N().inverse()

        if self.domain() is None or self.codomain() is None:
            raise ValueError(
                "Dual can only be computed once domain and codomain are known"
            )

        return LinearIsomorphism(domain = self.codomain(), codomain = self.domain(), N = N_inverse)

    def inverse_image(self, P):
        f=self.dual()
        return f(P)

    def apply_isomorphism(self, P):
        if self.N() is None:
            raise ValueError("Isomorphism matrix not specified")
        N = self.N()
        x, y, z, t = P.coords()
        X = N[0, 0] * x + N[0, 1] * y + N[0, 2] * z + N[0, 3] * t
        Y = N[1, 0] * x + N[1, 1] * y + N[1, 2] * z + N[1, 3] * t
        Z = N[2, 0] * x + N[2, 1] * y + N[2, 2] * z + N[2, 3] * t
        T = N[3, 0] * x + N[3, 1] * y + N[3, 2] * z + N[3, 3] * t

        return (X, Y, Z, T)

    def image(self, P):
        #print(self)
        return self.apply_isomorphism(P)

    def fold_compose(self, iso):
        assert iso.__class__ == Isomorphism
        domain = self.domain()
        codomain = iso.codomain()
        N1 = self.N()
        N2 = self.N()
        N = N2 * N1
        return Isomorphism(domain=domain, codomain=codomain, N=N)

class LinearChangeModel(LinearIsomorphism):
    pass

class LinearAutomorphism(LinearIsomorphism):
    def __init__(self, domain, **kwds):
        super().__init__(domain=domain, codomain=domain, **kwds)

class Automorphism(Isomorphism):
    def __init__(self, domain, **kwds):
        super().__init__(domain=domain, codomain=domain, **kwds)

class Translation(LinearAutomorphism):
    pass

class ChangeModel(Morphism):
    pass

class TrivialChangeModel(ChangeModel):
    @staticmethod
    def identity(K, image=None):
        if image is None:
            image = K
        elif image == "Id":
            def image(P):
                return P
        return TrivialChangeModel(domain=K, codomain=K, image=image, inverse_image=image)

    def __init__(self, domain = None, codomain = None, image = None, inverse_image = None, **kwds):
        if image is None:
            image = codomain
        if inverse_image is None:
            inverse_image = domain
        super().__init__(domain=domain, codomain=codomain, image=image, inverse_image=inverse_image,**kwds)

class Isogeny(Morphism):
    pass

class TranslatedIsogeny(Isogeny):
    pass
