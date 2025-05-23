AttachSpec("../Galois-Frob-Polys/src/spec");
AttachSpec("../src/spec");

// Records for Stab(w)-conjugacy classes of subgroups 
ClassInfo := recformat< label, 
                        group, 
                        anglerank, 
                        geomsimple, 
                        admissible, 
                        transitive, 
                        exceptional >;

YesNo := func<flag | flag select "Yes" else "No">;

function IsInteresting(G, wt)
  // A weighted permutation representation is "interesting" for the purposes
  // of this file if it is geometrically simple, strongly p-admissible (for some
  // p) and has non-maximal angle rank
  d := #wt div 2;
  ar := AngleRank(G, wt);
  
  if ar eq d then
    return false, ar;
  elif not IsGeometricallySimple(G, wt) then
    return false, ar;
  elif not HasAdmissibleFiltration(G, wt) then
    return false, ar;
  end if;
  return true, ar;
end function;

function GetInteresting(wt, GG, labs : magma_backup:=true)
  // Returns a list of all "interesting" (as defined above) Stab(w)-conjugacy
  // classes 
  d := #wt div 2;
  
  // housework
  disp_wt := [RealField(3)!w : w in wt];
  disp_wt := [Sprint(w) cat "_" : w in disp_wt[1..#wt-1]] cat [Sprint(disp_wt[#wt])];
  disp_wt := &cat disp_wt;
  backup_file := Sprintf("output-data/%o-%o.m", d, disp_wt);
  
  cclasses := {};
  all_subs := &cat GG;
  all_labs := &cat labs;
  
  for G in all_subs do
    Include(~cclasses, WeightConjugacyClass(PartitionFromSlopes(wt), G));
  end for;
  cclasses := [[G : G in cl] : cl in cclasses]; // All Stab(w)-ccs of subgroups
  
  cclasses_labs := [];
  for cl in cclasses do
    this_cl := [];
    for G in cl do
      i := Index(all_subs, G);
      Append(~this_cl, all_labs[i]);
    end for;
    Append(~cclasses_labs, this_cl);
  end for;

  ret := [];

  for cl in cclasses_labs do
    lab := cl[1];
    i := Index(all_labs, lab);
    G := all_subs[i];
    
    flag, ark := IsInteresting(G, wt);
    
    if flag then
      isexceptional := IsExceptional(G, wt);
      Append(~ret, rec<ClassInfo | 
                      label:=lab, 
                      group:=G, 
                      anglerank:=ark, 
                      geomsimple:=true, 
                      admissible:=true, 
                      transitive:=true, 
                      exceptional:=isexceptional>);
    end if;
    
    // Print to backup file
    if (Integers(100)!Index(cclasses_labs, cl) eq Integers(100)!#cclasses_labs) 
       and magma_backup then
      ret_info := [
        <r`label, 
          [Eltseq(g) : g in FewGenerators(r`group)], 
          r`anglerank, 
          r`exceptional> : 
        r in ret
      ];
      PrintFile(backup_file, "// label, group, angle rank, exceptional\n" : 
                Overwrite:=true);
      PrintFile(backup_file, ret_info);
    end if;
  end for;
  
  return ret;
end function;


function LoadInInfo(wt : directory:="output-data")
  // Load into a list the data stored in the output file generated by 
  // GetInteresting
  d := #wt div 2;
  W := W2(d);
  
  disp_wt := [RealField(3)!w : w in wt];
  disp_wt := [Sprint(w) cat "_" : w in disp_wt[1..#wt-1]] cat [Sprint(disp_wt[#wt])];
  disp_wt := &cat disp_wt;
  backup_file := Sprintf("%o/%o-%o.m", directory, d, disp_wt);
  info := eval Read(backup_file);
  
  ret := [];

  for eg in info do
    lab, gg, ark, isexceptional := Explode(eg);
    G := sub<W | [W!g : g in gg]>;
    Append(~ret, rec<ClassInfo | 
                    label:=lab, 
                    group:=G, 
                    anglerank:=ark, 
                    geomsimple:=true, 
                    admissible:=true, 
                    transitive:=true, 
                    exceptional:=isexceptional>);
  end for;
  
  return ret;
end function;

