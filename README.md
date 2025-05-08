This is the GitHub repository containing code related to the paper [*Galois groups of simple abelian varieties over finite fields and exceptional Tate classes*](https://arxiv.org/) by Santiago Arango-Piñeros, Sam Frengley, and Sameera Vemulapalli [arXiv:nnnnnn](https://arxiv.org/), [DOI:nnnnnn](https://www.doi.org/).

The code in this repository is written in [`Magma`](http://magma.maths.usyd.edu.au/magma/).

## Directory structure
- `src` contains the main intrinsics provided (in particular an implementation of Algorithm 6.2 and code to check in a weighted permutation representation is exceptional).

- `scripts` contains Magma scripts for verifying claims made in the article.

- `1-4-theorem-data` contains a list of exceptional weighted permutation representations in dimensions 2, 3, 4, 5, and 6.

- `Galois-Frob-Polys` (the [code](https://github.com/sarangop1728/Galois-Frob-Polys.git) associated to our previous joint article) is included as a submodule.