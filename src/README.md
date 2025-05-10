This directory contains three Magma files, each containing several intrinsics. Examples of the use of these intrinsics can be found in `../scripts/6-2-example.m` and `../scripts/6-2-algorithm-example.m`.

## `exceptionality.m`

The main intrinsics provided by this file are 
- `IsGeometricallySimple(G, w)` taking as input a subgroup $G \subset W_{2d}$ and a weight function $w$ and returning `true` if and only if $(w,G)$ is a geometrically simple weighted permutation representation
- `HasAdmissibleFiltration(G, w)` with input as above and returning `true` if and only if $(w,G)$ has a strongly $p$-admissible filtration for some prime number $p$.
- `IsExceptional(G, w)` with input as above and return `true` if and only if $(w, G)$ is exceptional.

## `numerical.m`

This file provides *non-rigorous* sanity checks. In particular, we employ the LLL algorithm on complex approximations to (normalised) Frobenius eigenvalues, in order to *numerically* compute the angle rank of a $q$-Weil polynomial and whether or not the corresponding abelian variety (assumed geometrically simple with commutative endomorphism algebra) carries exceptional Tate classes. The main intrinsics are:
- `NumericallyExceptionalRelations(P)` given a $q$-Weil polynomial $P$, determines numerical approximations to the (nontrivial) relations of the form $\prod \lambda_{\mathsf{i}} \prod \bar{\lambda}_{\mathsf{j}} = \zeta$
- `IsNumericallyExceptional(P)` input as above, returns `true` if $P$ gives rise to an exceptional weighted permutation representation (up to numerical tolerance), as per Lemma 4.3.

## `algo.m`

This file provides a reference implementation of Algorithm 6.2 from the paper. The main intrinsic is `wprToWeilPolynomial` and can be called in three ways (the second being the most common):

1. `wprToWeilPolynomial(w,K,p)` takes a weight $w$, a CM field $K$, and a prime number $p$ and returns a $p^k$-Weil polynomial provided that there exists a prime $\mathfrak{A}$ above $p$ in $L$ (the Galois closure of $K$) with strongly $p$-admissible ramification filtration for ***some $W_{2d}$-conjugate*** of $Gal(L/\mathbb{Q})$ in $W_{2d}$.

2. `wprToWeilPolynomial(G,w,K,p:EditG=false)` takes a subgroup $G \subset W_{2d}$, a weight $w$, a CM field $K$ with Galois group $W_{2d}$-conjugate to $G$, and a prime number $p$ and returns a $q = p^k$-Weil polynomial provided that there exists a prime $\mathfrak{A}$ above $p$ in $L$ (the Galois closure of $K$) with strongly $p$-admissible ramification filtration with respect to the weighted permutation representation $(w, G)$. If the flag `EditG` it set to `true` the group $G$ may be replaced by a $W_{2d}$-conjugate (possibly chaning the weighted permutation representation).

3. `wprToWeilPolynomial(G,w,K,L,P:EditG=false)` is the same as (2) except the user provides the Galois closure $L/K$ and the prime `P` ($=\mathfrak{A}$) in $\mathcal{O}_L$.

A verbose flag `Weil` is provided (taking values between $0$ to $3$). 

### Possible improvements

The reference implementation given here has three main limitations.
- It involves computing the *explicit* Galois action on a splitting field $L/K$. The degree of $L$ can be very large (even in dimension $5$) and this is the central limiting factor. It is possible that this can be avoided via either of the following:
    - Working instead with a $p$-adic approximation to this action (maybe via `Magma`'s `GaloisGroup` machinery when $p$ is unramified in $K$)
    - Making computations only in $K$.
- The implimentation involves computing the full ring of integers of $L$. It should also be possible to avoid this.
- Finally, a smaller issue is that the ideal intersection algorithm `ComputeIdealIntersection` involves some randomness (essentially it takes traces and norms of random ideal elements). Sometimes this can lead to dramatic slow downs.