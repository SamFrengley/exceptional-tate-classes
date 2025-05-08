// We have precomputed for you the output of `1-4-theorem-data-generate.m` (this
// is the content of the directory `output-data`. However, to believe us you may 
// want to run that script yourself.

// Reminder, `output-data` contains all weighted permutation representations (up
// to Stab(w)-conjugacy) such that they are
//    * geometrically simple
//    * strongly p-admissible for some p (actually, not quite
//       we only check the conditions in Remark 2.10)
//    * non-maximal angle rank
// it only suffices to check that the only non-exceptional wpr with angle rank 
// < 6 is this one!

load "data-generation.m";
load "newton-polys.m";

full_list := [];

for wt in wt6 do
  all_interesting := LoadInInfo(wt);
  for my_case in all_interesting do
    if (my_case`anglerank lt 6) and (not my_case`exceptional) then
      Append(~full_list, my_case);
    end if;
  end for;
end for;

assert #full_list eq 1;
full_list;
