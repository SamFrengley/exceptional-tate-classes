load "data-generation.m";
load "newton-polys.m";

needs_ramification := [];

for wts in all_wts do
  for wt in wts do
    table := LoadInInfo(wt);
    for row in table do
      G := row`group;
      flag, filtrs := HasAdmissibleFiltration(G, wt);
      assert flag;
      if not exists(fil){
                     fil : fil in filtrs | IsCyclic(fil[1]) and 
                                           #fil[2] eq 1
                   } then
         Append(~needs_ramification, <wt,row>);
         if forall{fil : fil in filtrs | #fil[3] ne 1} then
           ps := {};
           for fil in filtrs do
             _, p := IsPrimePower(#fil[3]);
             Include(~ps, p);
           end for;
           
           printf "With \n"
                  cat "- dimension %o,\n"
                  cat "- NP %o,\n"
                  cat "- WPR %o\n"
                  cat "we are in the worst case. Each admissible subgoup\n"
                  cat "need wild ramification. Possible primes are %o.\n\n",
                  (#wt div 2), wt, row`label, ps;
         end if;
      end if;
    end for;
  end for;
end for;
