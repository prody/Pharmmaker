# 2019-06-30 Bahar LAB, University of Pittsburgh  bahar@pitt.edu
# Written by Ji Young Lee
# Write PDB files for selected snapshots with dominant interactions from
# Druggability molecular dynamics simulations

set strucs APDB
set traj_list ADCD
set FRAME AAA
set PROBE BBB
set RESID CCC
set trajNum FF

# Take parameter values from input arguments as far as possible
for {set index 0} {$index < $argc -1} {incr index} {
  if {$index eq  0} {set strucs [lindex $argv $index]}
  if {$index eq  1} {set traj_list [lindex $argv $index]}
  if {$index eq  2} {set FRAME [lindex $argv $index]}
  if {$index eq  3} {set PROBE [lindex $argv $index]}
  if {$index eq  4} {set RESID [lindex $argv $index]}
  if {$index eq  5} {set trajNum [lindex $argv $index]}
}

set trajCount 0
foreach trajFile $traj_list {

  if { [llength $strucs] > 1 } {
    set struc [lindex $strucs $dcdNum]
  } else {
    set struc $strucs
  }

  if {$trajCount eq $trajNum} {
    mol load pdb $struc 
    mol addfile  $trajFile first 0 last -1 step 1 waitfor -1

    animate write pdb pro-$trajNum-$FRAME.pdb beg $FRAME end $FRAME waitfor all sel [atomselect top "protein"]
    animate write pdb lig-$trajNum-$FRAME-$PROBE-$RESID.pdb beg $FRAME end $FRAME waitfor all sel [atomselect top "resname $PROBE and resid $RESID"]
  }

  incr trajCount
}

exit
