declare verbose Weil, 3;

W2dToS2dLabel := func<i,d | Index(X2(d), i)>;
S2dToW2dLabel := func<i,d | X2(d)[i]>;


intrinsic WeightStabiliser(G::GrpPerm, wt::SeqEnum) -> GrpPerm
{The subgroup of G stabilising the weight function w}
  d := #wt div 2;

  // Implimentation quirk   
  if forall{i : i in [1..#wt-1] | wt[i] le wt[i+1]} then
    wt := wt[1..d] cat Reverse(wt[d+1..2*d]);
  end if;
  
  wt_partition := [];
  slopes := {@ w : w in wt @};
  for s in slopes do
    inds := {S2dToW2dLabel(i, d) : i in [1..2*d] | wt[i] eq s};
    Append(~wt_partition, inds);
  end for;
  return Stabilizer(G, wt_partition);
end intrinsic;


intrinsic IsStronglyPAdmissible(filtr::SeqEnum[GrpPerm], 
                        w::SeqEnum, 
                        p::RngIntElt) -> BoolElt
{Checks if the filtration can be strongly p-admissible. 
Input is a list filtr=[G, D, G_0, G_1], a weight function and a prime p.}
  d := #w div 2;
  
  for i in [1..3] do
    if not filtr[i+1] subset filtr[i] then
      return false;
    end if;
  end for;  
  
  X2d := GSet(filtr[1]);
  H := WeightStabiliser(filtr[1], w);
  
  if not filtr[2] subset H then
    return false;
  end if;
  
  S := Orbits(filtr[2], X2d);
  S := [Set(s) : s in S];
  sums := [ &+[w[W2dToS2dLabel(i, d)] : i in s] : s in S ];  
  if not forall{s : s in sums | Denominator(s) eq 1} then
    return false;
  end if;
  
  if not IsNormal(filtr[2], filtr[3]) then
    return false;
  end if;
  if not IsCyclic(filtr[2] / filtr[3]) then
    return false;
  end if;
  if not IsNormal(filtr[3], filtr[4]) then
    return false;
  end if;
  if not IsCyclic(filtr[3] / filtr[4]) then
    return false;
  end if;
  if (Index(filtr[3], filtr[4]) mod p) eq 0 then
    return false;
  end if;
  
  if #filtr[4] ne 1 and (#filtr[4] mod p) ne 0 then
    return false;
  end if;
  
  return true; 
end intrinsic;


intrinsic NewtonSlopes(f::RngUPolElt, q::RngIntElt) -> SeqEnum
{Computes the (normalised wrt q) Newton slopes of f}
  _, p, k := IsPrimePower(q);
  cc := Reverse(Coefficients(f));
  vals := [Valuation(c, p) / k : c in cc];
  
  ret := [];
  done := false;
  our_vals := vals;
  while not done do
    slopes := [our_vals[i]/(i-1) : i in [2..#our_vals]];
    slope_i := Min(slopes);
    len := Max([i : i in [1..#slopes] | slopes[i] eq slope_i]);
    ret cat:= [slope_i : i in [1..len]];
    our_vals := [v - our_vals[len+1] : v in our_vals[len+1..#our_vals]];
    done := our_vals eq [0];
  end while;
  
  return ret;
end intrinsic;


function ClassGroupOrder(I)
  /*
  Compute the order of I in the class group via the maximally naive approach. 
  Second return argument is a principal generator for I^k.
 */
  k := 1;
  while true do
    if IsPrincipal(I^k) then
      _, alpha := IsPrincipal(I^k);
      return k, alpha;
    end if;
    k +:= 1;
  end while;
end function;


function AllRootsOfUnity(K)
  /*
  Returns the set of roots of unity of K (not +- 1)
 */
  // Only allowed those with degree dividing degree of K. Degree of prim root of 
  // unity is phi(m) >= sqrt(m/2). Thus 2*deg(K)^2 >= m
  n := Degree(K);
  max_m := 2*n^2;
  allowed_m := [ m : m in [3..max_m] | n mod m eq 0 ];
  big_phi := &*[CyclotomicPolynomial(m) : m in allowed_m];
  return [r[1] : r in Roots(big_phi, K)];
end function;


function permutationFromOrderedRoots(g,tau,roots)
  /*
  Given an element of the Galois group G (which acts on the set of roots),
  return the permutation corresponding to the ordering of the roots.
  */
  return [Index(roots, tau(g)(r)): r in roots];
end function;


function conjugateToCorrectComplexConj(cc)
  /*
  Input: cc = (1,2)(3,5)(4,6) in S6.
  Output: Permutation s such that s*cc*s^(-1) = (1,6)(2,5)(3,4).
  */
  n := Degree(Parent(cc));
  Sn := SymmetricGroup(n);
  d := n div 2;
  s := Identity(Sn);
  for j in [1..d] do
    l := [1..n]; l[d+j] := j^cc; l[j^cc] := d+j;
    t := Sn!l;
    s := t*s;
    cc := t*cc*t^(-1);
  end for;
  return s;
end function;


function PermutationToW2d(perm)
  n := Degree(Parent(perm));
  d := n div 2;
  X := X2(d);
  W := W2(d);
  return W!{@ X[i^perm] : i in [1..n] @};
end function;


function getGaloisRamData(L, roots, P)
  n := #roots;
  d := n div 2;
  
  W := W2(d);

  // Automorphism Group
  G, _, tau := AutomorphismGroup(L); 

  OL := Order(P);
  D := DecompositionGroup(P);
  ram0 := RamificationGroup(P,0); // Inertia
  ram1 := RamificationGroup(P,1);
  filtr := [G,D,ram0,ram1];
  
  // Conjugate into W(2d)
  Sn := SymmetricGroup(n);
  gens := [];
  for g in Generators(G) do
    Append(~gens, <g,permutationFromOrderedRoots(g,tau,roots)>);
  end for;
  G1 := sub<Sn | [g[2] : g in gens]>;
  phi1 := hom<G -> G1 | gens>;
  phi1_inv := hom<G1 -> G | [<phi1(h), h> : h in G]>;
  filtr := [phi1(GG) : GG in filtr];
  
  _,iota := HasComplexConjugate(L);
  cc := G1![Index(roots, iota(r)) : r in roots]; 
  s := conjugateToCorrectComplexConj(cc);

  phi2 := hom<G1 -> Sn | h:->s*h*s^(-1)>;
  G2 := Image(phi2);
  phi2 := hom<G1 -> G2 | h:->s*h*s^(-1)>;
  phi2_inv := hom<G2 -> G1 | [<phi2(h), h> : h in G1]>;
  cc := phi2(cc); 
  filtr := [phi2(GG) : GG in filtr];
  
  roots := [roots[i^s]: i in [1..n]];
  
  // Coerce H to W2(d).
  W := W2(d);
  G3 := sub<W | [PermutationToW2d(h): h in Generators(G2)]>;
  phi3 := hom<G2 -> G3 | h:->PermutationToW2d(h)>; 
  phi3_inv := hom<G2 -> G3 | [<phi3(h), h> : h in G2]>; 
  filtr := [phi3(GG) : GG in filtr];
  cc := phi3(cc);
  
  // phi_inv := Inverse(phi1*phi2*phi3);
  phi_inv := phi3_inv*phi2_inv*phi1_inv;
  tau2 := map<G3 -> Codomain(tau) | g:->tau(phi_inv(g))>;
  
  return filtr, G, phi_inv, tau2;
end function;


function getConjugateGaloisRamData(G, L, roots, P)
  d := #roots div 2;
  W := W2(d);
  
  filtr, AutL, H_to_AutL, tau := getGaloisRamData(L, roots, P);
  H := filtr[1];
  
  _, s := IsConjugate(W, H, G); 
  phi := hom<H -> G | [s^(-1)*h*s : h in Generators(H)]>; 
  filtr := [phi(GG) : GG in filtr];
  
  phi_inv := Inverse(phi);
  tau2 := map<G -> Codomain(tau) | g:->tau(phi_inv(g))>;

  return filtr, AutL, phi_inv*H_to_AutL, tau2;
end function;


function ComputeIdealIntersection(I, P_divide_I, OK, H, tau)
  // I an ideal of L, K = L^H
  J := ideal<OK | [&*[tau(g)(i) : g in H] : i in Generators(I)] cat 
                  [&+[tau(g)(i) : g in H] : i in Generators(I)]>;
  done := false;

  while not done do
    Iprime := ideal<Order(I) | Generators(J)>;
    if forall{P : P in P_divide_I | Valuation(Iprime, P) eq Valuation(I, P)} then
      done := true;
    else
      elt := &+[Random(-10,10)*i : i in Generators(I)];
      a,b := TwoElement(J);
      J := ideal<OK | [a,b] cat 
                      [&*[tau(g)(elt) : g in H], &+[tau(g)(elt) : g in H]]>;
    end if;
  end while;
  return J;
end function;


function ComputeIdealNorm(I, OK, H, tau)
  // I an ideal of L, K = L^H
  J := ideal<OK | [&*[tau(g)(i) : g in H] : i in Generators(I)]>;
  return J;
end function;


function GetCorrectPrimeProduct(G_act_on_P, weight_exponents)
  // Input is a pair of lists of the same length. The first is the 
  // orbits g(P) for g \in G and the latter is their exponents w(g^(-1)(1))
  
  all_P := Set(G_act_on_P);
  
  our_product := [<P, weight_exponents[[i : i in [1..#G_act_on_P] | G_act_on_P[i] eq P]]> :
                  P in all_P];
  our_product := [<P[1], &+P[2]> : P in our_product];                  // the factors to be 
                                                                       //     multiplied with
                                                                       //     exponents
  c_num := LCM([Denominator(P[2]) : P in our_product]);                // numerator of c
  our_product := [<P[1], c_num*P[2]> : P in our_product];              // the factors to be 
                                                                       //     multiplied
                                                                       //     with c_num*exp
  c_den := GCD([Numerator(P[2]) : P in our_product]);                  // denominator of c
  our_product := [<P[1], Integers()!(P[2]/c_den)> : P in our_product]; // the factors to be 
                                                                       //     multiplied
                                                                       //     with exp/c_den  
  c := c_num/c_den;                                                    // c in the paper
  I := &*[P[1]^P[2] : P in our_product];
 
  return I, c;
end function;


function GetExponente(I, all_P, OK, H, tau)
  /*
  Gets the exponent e
  */ 
  // ramification index is multiplicative in towers.
  ram_L_over_K := {};
  for Q in all_P do
    q := Factorisation(ComputeIdealNorm(Q, OK, H, tau))[1][1];
    e_q := RamificationIndex(q);
    e_Q := RamificationIndex(Q);
    Include(~ram_L_over_K, <q, e_Q div e_q, Valuation(I, Q)>);
  end for;
  
  e_i := [ LCM(ram[2], ram[3]) div ram[3] : ram in ram_L_over_K | ram[3] ne 0];
  e := LCM(e_i);
  
  vprintf Weil, 2 : "CHECK: It is %o that e_q | e*v_q(I) for \n"
                    cat "       all q -- this should be true.\n\n", 
                    forall{ram : ram in ram_L_over_K | e*ram[3] mod ram[2] eq 0};
  return e;  
end function;


function FindAppropriateP(filtr, AutL, G_to_AutL, tau, P, w)
  _, p := IsPrimePower(Generators(P meet IntegerRing())[1]);
  
  Lcosets := {@ {@ g*h : h in filtr[2] @} : g in filtr[1] @};
  Lcosets := {@ coset[1] : coset in Lcosets @};                       // Reps
  
  for rep in Lcosets do
    new_P :=  tau(rep)(P);
    new_filtr := [Conjugate(H, rep) : H in filtr];
    
    // Debug check
    vprintf Weil, 3 : "CHECK: Conjugated to correct decomp group, %o.\n\n", 
                     forall{g : g in new_filtr[2] | tau(g)(new_P) eq new_P};
   
    
    if IsStronglyPAdmissible(new_filtr, w, p) then
      return new_filtr, new_P;
    end if;
  end for;

  error "No choice of prime above p gives p-admissible filtration";
end function;


intrinsic wprToWeilPolynomial(G::GrpPerm,
                              weight::SeqEnum,
                              K::FldNum,
                              L::FldNum,
                              P::RngOrdIdl : 
                              EditP:=false) -> RngUPolElt, RngIntElt
{On input a WPR group G and a CM field K, splitting field L (s.t Gal(L) = G), a
weighting w, and a prime ideal P of L, return a Weil polynomial which defines K
as a number field. If the flag EditP is true then algorithm will run through all
divisors P of p in O_L and check p-admissibility.}
  
  d := Degree(K) div 2;
  w := weight[1..d] cat Reverse(weight[d+1..2*d]); // Implementation quirk. To see why look
                                                   //    at the ordering of elts of X2(d).
  _, p := IsPrimePower(Generators(P meet IntegerRing())[1]);
  OL := MaximalOrder(L);

  roots := [r[1] : r in Roots(DefiningPolynomial(K), L)];
  filtr, AutL, G_to_AutL, tau := getConjugateGaloisRamData(G, L, roots, P);
  
  vprintf Weil, 1: "INFO: The Galois group is the group isomorphic to %o and\n"
                   cat "      the filtration is %o\n\n", 
                   GroupName(G), 
                   [GroupName(fil) : fil in filtr];
  
  if not EditP then
    require IsStronglyPAdmissible(filtr, w, p) : "This is not a strongly p-admissible filtration";
  else 
    filtr, P := FindAppropriateP(filtr, AutL, G_to_AutL, tau, P, w);
  end if;

  /* [tau(g)(P) eq P : g in filtr[2]]; */

  H := Stabilizer(G, 1); 
  H2 := [G_to_AutL(h) : h in Generators(H)];
  H2 := sub<AutL | H2>;
  new_K := FixedField(L, H2);
  OK := Integers(new_K);
  _, my_phi := IsSubfield(new_K, L);
  
  gg := [g : g in G];
  weight_exponents := [w[ W2dToS2dLabel(1^(g^(-1)), d) ] : g in gg];
  G_act_on_P := [ tau(g)(P) : g in gg ];
  I, c := GetCorrectPrimeProduct(G_act_on_P, weight_exponents);
  e := GetExponente(I, Set(G_act_on_P), OK, H, tau);
  J := ComputeIdealIntersection(I^e, Set(G_act_on_P), OK, H, tau);
  
  //----------------- Verbose print of sanity checks --------------------------//
  vprintf Weil, 2 : "CHECK: It is %o that L^H = K -- this should be true.\n\n",
                    IsIsomorphic(K, new_K);
                                
  vprintf Weil, 2 : "CHECK: It is %o that I^H = H -- this should be true.\n\n", 
          forall{h : h in H | tau(h)(I) eq I};
  
  vprintf Weil, 1: "CHECK: It is %o that v(g(I)) = #D*c*w(g(1)) for all\n"
                   cat "       g in G -- this should be true.\n\n",
                   forall{g : g in G | 
                          Valuation(tau(g)(I), P) eq #filtr[2]*c*w[W2dToS2dLabel(1^(g), d)]
                         };
  
  JOL := ideal<OL | [j : j in Generators(J)]>;
  vprintf Weil, 2 : "CHECK: It is %o that J = OK \cap I.\n\n", 
                   forall{gP : gP in Set(G_act_on_P) | 
                          Valuation(I, gP) eq Valuation(JOL, gP)
                         };
  //---------------------------------------------------------------------------//
  
  n, alpha := ClassGroupOrder(J);
  
  _, cc := HasComplexConjugate(new_K);

  Nm_al := Integers()!Norm(alpha); 
  Nm_p := Integers()!Norm(new_K!p); 
  k := 2*Valuation(Nm_al,p) div Valuation(Nm_p,p); 
  u := alpha*cc(alpha) / p^k;
  
  // Compute maximal totally real subfield
  Kplus := Subfields(new_K, d);
  assert exists(Kplus){KK : KK in Kplus | IsTotallyReal(KK[1])};
  Kplus, emb := Explode(Kplus);
  
  if IsSquare(Kplus!u) then
    ell := 1;
    _, v := IsSquare(Kplus!u);
  else
    ell := 2;
    v := u;
  end if;
  
  pi := alpha^ell/v;
  
  // Now the final iteration
  for zeta in [1] cat AllRootsOfUnity(new_K) do
    f := MinimalPolynomial(zeta*pi);
    if Degree(f) eq 2*d then
      return PolynomialRing(Integers())!f, 
             p^(Valuation(TrailingCoefficient(f), p) div d);
    end if;
  end for;
  
  error "Algorithm failed to terminate with something good";
end intrinsic;


intrinsic wprToWeilPolynomial(G::GrpPerm,
                              w::SeqEnum,
                              K::FldNum,
                              p::RngIntElt:
                              EditG:=false) -> RngUPolElt, RngIntElt
{On input a WPR group G and a CM field K which Gal(L) = G, a weighting w, and a
prime number p, return a Weil polynomial which defines K as a number field. If
EditG is true, then the algorithm treats G only as a conjugacy class of subgroups
of W_2d, rather than a w-conjugacy class.}
  
  d := Degree(K) div 2;
  L, roots := SplittingField(K);
  OL := MaximalOrder(L);
  
  P := Factorisation(ideal<OL | p>)[1][1];
  
  if not EditG then
    return wprToWeilPolynomial(G, w, K, L, P : EditP:=true);
  else
    W2d := W2(d);
    class_G := Class(W2d, G);
    for new_G in class_G do
      try 
        return wprToWeilPolynomial(G, w, K, L, P : EditP:=true);
      catch e
        assert 0 eq 0;
      end try;
    end for;
    error "None of the conjugates of G work";
  end if;
end intrinsic;


intrinsic wprToWeilPolynomial(w::SeqEnum,
                              K::FldNum,
                              p::RngIntElt) -> RngUPolElt, RngIntElt
{On input a CM field K which Gal(L) = G, a weighting w, and a prime number p,
return a Weil polynomial which defines K as a number field.}
  G := GaloisGroup(K);
  d := Degree(K) div 2;
  W := W2(d);
  gps := [H`subgroup : H in Subgroups(W : OrderEqual:=#G) | 
          IsIsomorphic(H`subgroup, G) and 
          IsTransitive(H`subgroup)];
  gps := &join { Class(W, H) : H in gps };

  for new_G in gps do
    try 
      return wprToWeilPolynomial(new_G, w, K, p);
    catch e
      assert 0 eq 0;
    end try;
  end for;
  error "None of the subgroups H work";
end intrinsic;
