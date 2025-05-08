declare verbose Weil, 3;

intrinsic NumericalRelations(P::RngUPolElt) -> Any
{}
  d := Degree(P) div 2;
  CC := ComplexField(1200);
  PP<x> := PolynomialRing(CC);
  P := PP!P;
  
  rr := [r[1] : r in Roots(P)];
  rr := [r : r in rr];
  rs := [];
  for r in rr do
    if not exists{x : x in rs | Abs(Re(x)-Re(r)) lt 10^(-20)} then
      Append(~rs, r);
    end if;
  end for;

  rs := [Arg(r) / (2*Pi(CC)) : r in rs];
  
  // base change torsion with large number
  tors_kill := &*[p^2 : p in PrimesInInterval(1, 2*d)];
  rs := [tors_kill*r : r in rs];
  
  rs := rs cat [1]; // because argument of 2*pi is the same as 0
  
  B := Basis(AllLinearRelations(rs, 400));
  B := [Eltseq(b)[1..d] : b in B];
  return B;
end intrinsic;

intrinsic NumericalExceptionalRelations(P::RngUPolElt) -> SeqEnum, RngIntElt
{Returns all the excpetional relations and the minimal codimension, r}
  d := Degree(P) div 2;
  L := NumericalRelations(P);
  L := Lattice(Matrix(Integers(), #L, #L[1], L));

  vv := ShortVectors(L, d); // max length is length of (1,...,1) = sqrt(d) 
  vv := [Eltseq(v[1]) : v in vv];

  // sign relations
  vv := [v : v in vv | (Set(v) meet {-1,0,1}) eq Set(v)];
  
  // even number of pm 1
  vv := [v : v in vv | #[x : x in v | x in {-1,1}] mod 2 eq 0];
  
  if #vv eq 0 then
    r := -1;
  else
    r := Min([#[x : x in v | x in {-1,1}] div 2 : v in vv]);
  end if;
  
  return vv, r;
end intrinsic;


intrinsic IsNumericallyExceptional(P::RngUPolElt) -> BoolElt
{Checks the condition for being exceptional}
  return #NumericalExceptionalRelations(P) ne 0;
end intrinsic;
