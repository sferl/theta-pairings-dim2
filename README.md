# Hyperelliptic biextension pairings

Author: Alessandro Sferlazza. <br>
Repository largely based on:
- [An Algorithmic Approach to (2, 2)-isogenies in the Theta Model](https://github.com/ThetaIsogenies/two-isogenies) (Pierrick Dartois, Luciano Maino, Giacomo Pope, and Damien Robert)
- [Kummer Line](https://gitlab.inria.fr/roberdam/kummer-line) (Damien Robert)

This repository wants to be a SageMath implementation of pairings on 2-dimensional abelian varieties (e.g., hyperelliptic Jacobians) using theta coordinates.

The algorithms are based on [Fast pairings via Biextensions](https://eprint.iacr.org/2024/517) (Damien Robert, 2024).

The code constitutes an early-stage proof-of-concept implementation, in need of further work.
My personal contribution is minimal and consists in having merged the two aforementioned repositories:
- [(2, 2)-isogenies](https://github.com/ThetaIsogenies/two-isogenies) provides an interface for 2-dimensional abelian varieties in the theta model. Alongside the existing projective arithmetic, we implemented the affine *cubical* arithmetic described in [Fast pairings via Biextensions](https://eprint.iacr.org/2024/517), less efficient than the former, but necessary for our pairing algorithms. Its code is contained in folders <code>theta_isogenies</code>, <code>theta_structures</code>, <code>utilities</code>
- [Kummer Line](https://gitlab.inria.fr/roberdam/kummer-line) implements [Fast pairings via Biextensions](https://eprint.iacr.org/2024/517), mainly via the <code>Biextension</code> class. Originally implemented to work with biextensions of elliptic curves (also in the theta model), the code also works for higher-dimensional abelian varieties. Code contained in folder <code>biextensions</code>.

Requirements: SageMath 10.3



