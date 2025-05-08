// We have precomputed for you the output of `1-4-theorem-data-generate.m` (this
// is the content of the directory `output-data`. However, to believe us you may 
// want to run that script yourself.

// Reminder, `output-data` contains all weighted permutation representations (up
// to Stab(w)-conjugacy) such that they are
//    * geometrically simple
//    * strongly p-admissible for some p (actually, not quite
//       we only check the conditions in Remark 2.10)
//    * non-maximal angle rank
// thus it only suffices to check that the exceptional ones are listed in the 
// directory `../1-4-theorem-data`

load "data-generation.m";
load "newton-polys.m";

for wts in all_wts do
  for wt in wts do
    all_interesting := LoadInInfo(wt);
    thm_data := LoadInInfo(wt : directory:="../1-4-theorem-data");
    for my_case in all_interesting do
      if my_case`exceptional then
        // we just need to show it is listed in the data
        assert my_case`label in [eg`label : eg in thm_data];
      end if;
    end for;
  end for;
end for;
