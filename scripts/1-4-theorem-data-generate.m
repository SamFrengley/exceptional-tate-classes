// This file computes all the nessicary data for every (non-supersingular) 
// Newton polygon in dimension <= 6. The nontrivial content here is in the file
// `data-generation.m`.

SetColumns(1000); // print width

load "data-generation.m";
load "newton-polys.m";

for wts in all_wts do
  d := #wts[1] div 2;
  GG, labs := AdmissibleSubgroupsAndLabels(d); // get all subgroups of W_2d 
                                               // (admissible is not the same 
                                               //  adjective as p-admissible, sorry)
  for wt in wts do
    _ := GetInteresting(wt, GG, labs); // outputs to a file in `output-data`
  end for;
end for;
