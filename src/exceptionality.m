declare verbose Weil, 3;

W2dToS2dLabel := func<i,d | Index(X2(d), i)>;
S2dToW2dLabel := func<i,d | X2(d)[i]>;

function ReorderWeight(wt)
  d := #wt div 2;
  if forall{i : i in [1..#wt-1] | wt[i] le wt[i+1]} then
    wt := wt[1..d] cat Reverse(wt[d+1..2*d]);
  end if;
  return wt;
end function;


intrinsic PartitionFromSlopes(wt::SeqEnum) -> SetEnum
{Takes Newton slopes and returns partitions}
  d := #wt div 2;
  wt := ReorderWeight(wt);
  slopes := {s : s in wt};
  ii := [ {i : i in [1..#wt] | wt[i] eq s } : s in slopes ];
  return [ {S2dToW2dLabel(i, d) : i in si} : si in ii ];
end intrinsic;


intrinsic IsGeometricallySimple(G::GrpPerm, wt::SeqEnum) -> BoolElt
{Checks if you are geometrically simple}
  
  d := #wt div 2;
  require 2*d eq #GSet(G) : "Group doesn't act on the right set?";
  
  // Implimentation quirk   
  wt := ReorderWeight(wt);

  phi1 := { <wt[ W2dToS2dLabel( 1^(g^(-1)), d) ], g> : g in G };
  
  Stab_phi1 := {g : g in G | {<s[1], s[2]*g> : s in phi1} eq phi1};
  
  return Stab_phi1 eq {s : s in Stabilizer(G, 1)};
end intrinsic;

function GenerationCondition(D, G0, G1) 
  H, q := D/G1;
  G0bar := q(G0);
  
  taus := {tau : tau in H | sub<H|tau> eq G0bar};
  sigtaus := &join{{[sigma,tau] : tau in taus | sub<H | [sigma,tau]> eq H} : sigma in H};

  if #G1 eq 1 then
    // if G1 trivial then p can be anything (coprime to #(G0/G1) )
    // thus only need sig*tau*sig^(-1) in the subgroup <tau>
    sigtaus := {st : st in sigtaus | st[1]*st[2]*st[1]^(-1) in sub<H|st[2]>};
    return #sigtaus ge 1;
  end if;
  
  _, p := IsPrimePower(#G1);
  sigtaus := {st : st in sigtaus | st[1]*st[2]*st[1]^(-1) eq st[2]^p};

  return #sigtaus ge 1;  
end function;

intrinsic HasAdmissibleFiltration(G::GrpPerm, w::SeqEnum) -> BoolElt, SeqEnum
{Returns true if and only if there exists an admissible filtration of G.}
  d := #w div 2;
  w := ReorderWeight(w);
  
  Wstab := WeightStabiliser(G, w);
  DD := Subgroups(Wstab);
  DD := [D`subgroup : D in DD];
  
  all := [];
  
  for D in DD do
    X2d := GSet(G);
    S := Orbits(D, X2d);
    S := [Set(s) : s in S];
    sums := [ &+[w[W2dToS2dLabel(i, d)] : i in s] : s in S ];  
    if forall{s : s in sums | Denominator(s) eq 1} then
      
      G0s := NormalSubgroups(D);
      G0s := [G0`subgroup : G0 in G0s];
      G0s := [G0 : G0 in G0s | IsCyclic(D / G0)];
      for G0 in G0s do
        G1s := NormalSubgroups(G0);
        G1s := [G1`subgroup : G1 in G1s];
        G1s := [G1 : G1 in G1s | IsNormal(D, G1)];
        G1s := [G1 : G1 in G1s | #G1 eq 1 or (IsPrimePower(#G1))];
        for G1 in G1s do
          if IsCyclic(G0 / G1) then
            if #G1 eq 1 then
              Append(~all, [D,G0,G1]);
            else
              _, p := IsPrimePower(#G1);
              if (Index(G0, G1) mod p) ne 0 then
                if GenerationCondition(D,G0,G1) then
                  Append(~all, [D,G0,G1]);
                end if;
              end if;
            end if;
          end if;
        end for;
      end for;
      
    end if;
  end for;
  
  return #all ne 0, all;
end intrinsic;


function ExceptionalCondition(g, w, S_plus, S_minus)
  d := #w div 2;

  w := ReorderWeight(w);  
  
  RHS := &+([0] cat [w[ W2dToS2dLabel( i^(g^(-1)), d) ] : i in S_plus]) 
           - &+([0] cat [ w[ W2dToS2dLabel( i^(g^(-1)), d) ] : i in S_minus]);
  
  LHS := (#S_plus - #S_minus) div 2;
  return RHS eq LHS;
end function;


intrinsic IsExceptional(G::GrpPerm, w::SeqEnum) -> BoolElt, SeqEnum
{Returns true iff (G,w) is exotic wpr}
  d := #w div 2;

  w := ReorderWeight(w);

  S := {i : i in [1..d]};
  
  for S_plus in Subsets(S) do
    S_minuses := S diff S_plus;
    S_minuses := {S_minus : S_minus in Subsets(S_minuses) | 
                  (#S_minus + #S_plus) gt 0 and
                  (#S_minus + #S_plus) mod 2 eq 0 };

    for S_minus in S_minuses do
      if forall{g : g in G | ExceptionalCondition(g, w, S_plus, S_minus)} then
        return true, [S_plus, S_minus];
      end if;
    end for;
  end for;
  
  return false;
end intrinsic;
