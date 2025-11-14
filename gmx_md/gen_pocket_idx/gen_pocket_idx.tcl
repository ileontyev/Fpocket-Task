##  INPUT PARAMS
set pocketDir "../../fpocket-3_prot-atoms/7WSL_protprep_min_out/pockets"
set target_struct "pd1_7wsl_md298.pdb"
set nResShift 30

mol delete all
############################################################################
##               FORMATTED PRINTING
############################################################################
proc print_formatted {fp list maxlen} {
  set line {}
  set printed 0
  foreach word $list {
    set line1 "$line $word"
    set length [string length $line1]
    #puts "$line1 len=$length"
    if {$length <= $maxlen} {
      set line $line1
      set printed 0
    } else {
      puts $fp $line
      set printed 1
      set line " $word"
    }
  }
  # Print the last line if it hasn't been printed
  if {$printed != 1} {   
      puts $fp $line
  }
}

proc highlighting { colorId representation id selection } {
   puts "highlighting $id"
   mol representation $representation
   mol material "Diffuse" 
   mol color $colorId
   mol selection $selection
   mol addrep $id
}

proc highlighting { colorId representation id selection } {
   puts "highlighting $id"
   mol representation $representation
   mol material "Diffuse" 
   mol color $colorId
   mol selection $selection
   mol addrep $id
}
#######################################################################
#                         START
#######################################################################
set pocket_files [glob $pocketDir/pocket*_atm.pdb]
puts $pocket_files

## LOAD TARGET STRUC
set target_id [mol new $target_struct]
    mol delrep top $target_id
    highlighting Name "Lines" $target_id "protein and noh"
    highlighting Element "NewCartoon" $target_id "protein"
    #highlighting Name "Licorice" $target_id "not protein and not water"
    #highlighting "ColorID 7" "VdW 0.4" $target_id "protein and occupancy>0.95"
    display rendermode GLSL

set fp [open idx_pockets.ndx [list TRUNC CREAT WRONLY APPEND]]
set id 0
foreach pdbfile $pocket_files {
   ## READ POCKET ATOMS
   set str_id [mol new $pdbfile]
   puts "\n$pdbfile str_id=$str_id"
   set sel_pocket [atomselect top "all"]
   set nAtPocket [$sel_pocket num]]

   
   ## ATOM SELECTION by resid and atom names
   set resid_list [lsort -unique -integer [$sel_pocket get resid]]
   puts "residues: $resid_list"
   set sel_list {}
   foreach resid $resid_list {
      lappend sel_list "resid [expr $resid - $nResShift] and name [join [[atomselect top "resid $resid"] get name]]"
   }
   set sel_str [join $sel_list ") or ("]
   puts "sel_str: $sel_str"
   set sel_target_pocket [atomselect $target_id "($sel_str)"]
   set nAtPocket_target [$sel_target_pocket num]]
   if { $nAtPocket_target != $nAtPocket } {
       puts "ERROR: nAtPocket_target($nAtPocket_target) differs from nAtPocket ($nAtPocket)"
       return
   }
   puts "nAtPocket_target = $nAtPocket_target"
   
   ## WRITE INDICIES
   set idx_pocket_target [$sel_target_pocket get serial]
   set id [expr $id + 1]
   puts $fp "\[ pocket$id \]"
   print_formatted $fp $idx_pocket_target 80
   
   mol delete $str_id
   
   ## VISUALIZE POCKET ATOMS
   # mol material "Glass3"
   #highlighting "ColorID $id" "QuickSurf 1.0" $target_id "($sel_str)"
   highlighting "ColorID $id" "VDW 1.0" $target_id "serial $idx_pocket_target"
}
close $fp
display resetview
display rendermode GLSL
