This directory contains scripts which either give examples of the functionality of our code, or verify computational claims made in the article.

## Verifying claims

- The script `1-4-theorem-proof.m` proves Theorem 1.4 (under the assumption that `1-4-theorem-data-generate.m` has populated the data in `./output-data/`). 

- The script `1-4-theorem-data-generate.m` generates the data used in the proof of Theorem 1.4. This file takes ***a long time*** (several hours on a single core of one of the authors' laptops). For convinience we have saved this output in `./output-data/` (running the script will overwrite the saved data, though the output should be identical). The specific computations are:
    - We enumerate all (transitive) subgroups $G \subset W_{2d}$ up to $\mathrm{Stab}(w)$-conjugacy. 
    - We compute in each case for $\rho = (w,G)$:
        1. Whether $\rho$ is geometrically simple,
        2. The angle rank of $\rho$,
        3. Whether $\rho$ is *not* strongly $p$-admissible for any $p$ (using the criteria in Remark 2.10).

- The script `1-4-theorem-non-tame-remark.m` verifies the claim made in the paragraph directly following the statement of Theorem 1.4. Namely we check that there exists *exactly one* weighted permutation representation in dimension $\leq 6$ such that:
    - $\rho$ is geometrically simple
    - $\rho$ has angle rank $\delta_\rho < g$, 
    - $\rho$ is strongly $p$-admissible for some $p$, and 
    - every strongly $p$-admissible filtration of $\rho$ is wildly ramified.

- The script `6-6-remark-proof.m` verifies the claim that the weighted permutation representation in Section 6.2 (Corollary 1.5) is the only strongly $p$-admissible weighted permutation representation associated to a geometrically simple abelian $g$-fold with 
    - $g \leq 6$,
    - $\delta_A < g - 1$ (and $\delta_A < g$ if $g$ even), and
    - no exceptional Tate classes.

## Examples of the algorithm

The following code block generates the examples in Section 6.2.
```
AttachSpec("../Galois-Frob-Polys/src/spec");
AttachSpec("../src/spec");

// Define the number field
QQ := Rationals();
R<x> := PolynomialRing(QQ);
f := x^12 - 6*x^11 + 15*x^10 - 18*x^9 - 3*x^8 + 24*x^7 + 3*x^6 + 6*x^5 + 33*x^4 +
     152*x^3 + 273*x^2 + 198*x + 51;
K := NumberField(f);

// Define the weighting
w := [ 0, 0, 0, 1/3, 1/3, 1/3, 2/3, 2/3, 2/3, 1, 1, 1 ];

// Define the Galois group (as a subgroup of W_6)
W6 := W2(6);
G := sub<W6 | [ 
          W6!(2,4)(3,-5)(5,-3)(-2,-4), 
          W6!(1,-3,-4)(2,-6,5)(3,4,-1)(6,-5,-2)
        ]>;

// Get the examples
P2, q2 := wprToWeilPolynomial(G, w, K, 2);
P3, q3 := wprToWeilPolynomial(G, w, K, 3);
```
For more detailed examples:

- The script `6-2-example.m` illustrates the example in Section 6.2 (also Corollary 1.5). In particular it utilises some of the numerical algorithms in `../src/numerical.m` to give sanity checks on Corollary 1.5. 
    
- The script `6-2-algorithm-examples.m` provides several examples of how to use our implimentation of Algorithm 6.2. 


## Helper functions

- The file `data-generation.m` provides some helper functions for `1-4-theorem-data-generate.m` and functionality to read the stored data in `./output-data/` and `../1-4-theorem-data/` into `Magma`.

- The file `newton-polys.m` hardcodes the Newton polygons in dimensions $\leq 6$.
