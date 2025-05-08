AttachSpec("../Galois-Frob-Polys/src/spec");
AttachSpec("../src/spec");
SetVerbose("Weil", 0);

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


printf "The angle rank of rho=(w,G) is %o\n", AngleRank(G, w);
printf "It is %o that the wpr rho is geometrically simple\n",
        IsGeometricallySimple(G, w);
printf "It is %o that the wpr rho has a p-admissible filtration for some p\n\n", 
        HasAdmissibleFiltration(G, w);

/*
// Uncomment this block to generate for yourself the polynomials
// HINT: If it runs slowly, run it again because there is some randomness which
//       may slow you up.

P2, q2 := wprToWeilPolynomial(G, w, K, 2);
P3, q3 := wprToWeilPolynomial(G, w, K, 3);

printf "We found the %o-Weil polynomial\n   %o\n\n", q2, R!P2;
printf "We found the %o-Weil polynomial\n   %o\n\n", q3, R!P3;
*/

P2 := x^12 - 12*x^11 + 75*x^10 - 351*x^9 + 1392*x^8 - 4692*x^7 + 13912*x^6 
      - 37536*x^5 + 89088*x^4 - 179712*x^3 + 307200*x^2 - 393216*x + 262144;
q2 := 8;

P3 := x^12 - 3*x^11 + 14*x^9 - 21*x^8 - 27*x^7 + 120*x^6 - 81*x^5 -
       189*x^4 + 378*x^3 - 729*x + 729;
q3 := 3;

// Verify weighted permutation representation, angle rank, and exceptional
_, G := NewtonGaloisGroup(P2);
wt := NewtonSlopes(P2, q2);

printf "Numerically check that it is %o that P_2 is exceptional\n", 
       IsNumericallyExceptional(P2);
printf "Check that it is %o that P_2 is exceptional\n\n",
       IsExceptional(G, wt);

// Verify weighted permutation representation, angle rank, and exceptional
_, G := NewtonGaloisGroup(P3);
wt := NewtonSlopes(P3, q3);

printf "Numerically check that it is %o that P_3 is exceptional\n", 
       IsNumericallyExceptional(P3);
printf "Check that it is %o that P_3 is exceptional\n\n",
       IsExceptional(G, wt);
