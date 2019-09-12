# 2019-06-30 Bahar LAB, University of Pittsburgh  bahar@pitt.edu
# Written by Ji Young Lee
# Calculate interaction statistics from snapshots of 
# Druggability molecular dynamics simulations
# and rank dominant interactions


# Check cms
set struc      APDB
set dcd_list { ADCD }
set interval SSTEP 
set CUTOFF1 CUTOFF
set PROBE AAA
set RESID BBB
set CHAIN CCC
set TRJNUM RRR

# Take parameter values from input arguments as far as possible
for {set index 0} {$index < $argc -1} {incr index} {
  if {$index eq  0} {set struc [split [lindex $argv $index] ,]}
  if {$index eq  1} {set dcd_list [split [lindex $argv $index] ,]}
  if {$index eq  2} {set interval [lindex $argv $index]}
  if {$index eq  3} {set CUTOFF1 [lindex $argv $index]}
  if {$index eq  4} {set PROBE [lindex $argv $index]}
  if {$index eq  5} {set RESID [lindex $argv $index]}
  if {$index eq  6} {set CHAIN [lindex $argv $index]}
  if {$index eq  7} {set TRJNUM [split [lindex $argv $index] ,]}
}

set ofile [open ligbo-$TRJNUM.dat w]
set count 0

#write a header for the datafile 
#  puts $ofile "#frame    Calpha_dist    N-O_dist    N-O_dist    N-O_dist    N-O_dist    N-O_dist    N-O_dist"

foreach dcd_in $dcd_list {

  mol load pdb $struc
  mol addfile $dcd_in type dcd first 0 last -1 step $interval waitfor -1

  ##### alignment
  set ref_molid [molinfo top get id]
  set traj_molid [molinfo top get id]

  set calpha "alpha carbon"
  set CAreference [atomselect $ref_molid $calpha frame 0]
  set CAcompare [atomselect $traj_molid $calpha]
  set allsys "all"
  set allsysC [atomselect $traj_molid $allsys]
  ##### alignment


  set first_frame 0
  set num_frames [molinfo top get numframes] 

  for {set f $first_frame} {$f < $num_frames} {incr f} {
    incr count
    set time [expr ($count-1)*$interval]
    molinfo top set frame $f

    ##### alignment
    $CAcompare frame $f
    $allsysC frame $f
    set trans_matCA [measure fit $CAcompare $CAreference]
    $allsysC move $trans_matCA
    ##### alignment

    set gluOxy [atomselect top "protein and chain $CHAIN and resid $RESID and not hydrogen"]
    $gluOxy frame $f
    set Olist [$gluOxy list]

    ############ Check distance cutoff
    set argNit [atomselect top "resname $PROBE and not hydrogen"]
    $argNit frame $f
    set Nlist [$argNit list]
    set NNNN 0
    foreach atom1 $Olist {
      foreach atom2 $Nlist {
        set NOdist [measure bond [list [list $atom1] [list $atom2]]]
        if {$NOdist < $CUTOFF1} {
          set NNNN [expr $NNNN+1]
          set aato1 [expr $atom1+1]
          set aato2 [expr $atom2+1]
          set phe1 [atomselect [molinfo top get id] "serial $aato2"]
          set pha1 [atomselect [molinfo top get id] "serial $aato1"]
          set phe_vec1 [lindex [$phe1 get {resname}] 0]
          set phe_vec2 [lindex [$phe1 get {resid}] 0]
          set phe_vec3 [lindex [$phe1 get {name}] 0]
          set pha_vec1 [lindex [$pha1 get {resname}] 0]
          set pha_vec2 [lindex [$pha1 get {resid}] 0]
          set pha_vec3 [lindex [$pha1 get {name}] 0]
          animate write pdb v-$TRJNUM-$f.pdb beg $f end $f sel [atomselect top "protein or (resname $PROBE and resid $phe_vec2)"]
          puts $ofile "$f $aato1 $pha_vec1 $pha_vec2 $pha_vec3  $aato2 $phe_vec1 $phe_vec2 $phe_vec3  $NOdist"
        }
      }
    }	
    ############
  }
  #end loop over frames
  animate delete all
}
#end dtr list
flush $ofile
close $ofile
exit
