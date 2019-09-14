# 2019-06-30 Bahar LAB, University of Pittsburgh  bahar@pitt.edu
# Written by Ji Young Lee
# Write PDB files for selected snapshots with dominant interactions from
# Druggability molecular dynamics simulations

set strucs APDB
set traj_list ADCD
set interval SSTEP
set FRAME AAA
set PROBE BBB
set RESID CCC
set RRR FF

# Take parameter values from input arguments as far as possible
for {set index 0} {$index < $argc -1} {incr index} {
  if {$index eq  0} {set strucs [lindex $argv $index]}
  if {$index eq  1} {set traj_list [lindex $argv $index]}
  if {$index eq  2} {set interval [lindex $argv $index]}
  if {$index eq  3} {set FRAME [lindex $argv $index]}
  if {$index eq  4} {set PROBE [lindex $argv $index]}
  if {$index eq  5} {set RESID [lindex $argv $index]}
  if {$index eq  6} {set RRR [lindex $argv $index]}
}

foreach trajFile $traj_list {

  if { [llength $strucs] > 1 } {
    set struc [lindex $strucs $dcdNum]
  } else {
    set struc $strucs
  }

mol load pdb $struc 
mol addfile  $trajFile first 0 last -1 step $interval waitfor -1

animate write pdb pro-$RRR-$FRAME.pdb beg $FRAME end $FRAME waitfor all sel [atomselect top "protein"]
animate write pdb lig-$RRR-$FRAME-$PROBE-$RESID.pdb beg $FRAME end $FRAME waitfor all sel [atomselect top "resname $PROBE and resid $RESID"]

exit
