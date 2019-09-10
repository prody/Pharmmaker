# 2019-06-30 Bahar LAB, University of Pittsburgh  bahar@pitt.edu
# Written by Ji Young Lee
# Calculate high affinity residues with trajectory of 
# Druggability molecular dynamics simulations

namespace eval ::highAff:: {
  namespace export highAff

  ############ Check and adjust #################
  # PDB file name and path
  variable struc      APDB
  # DCD file name and path
  variable dcd_list { ADCD }
  # protein chain ID in pdb
  variable CHAIN      ACHAIN
  # residue number of protein for analysis
  variable RESID      ARESID
  # Probe molecule ID in pdb
  variable PROBE      APROBE
  # cutoff between residue and protein    
  variable CUTOFF     ACUTOFF
  # frequency of dcd for analysis
  variable STEP       ASTEP
  ############ Check and adjust #################
}

# Take parameter values from input arguments as far as possible
for {set index 0} {$index < $argc -1} {incr index} {
  if {$index eq  0} {set ::highAff::struc [split [lindex $argv $index] ,]}
  if {$index eq  1} {set ::highAff::dcd_in [split [lindex $argv $index] ,]}
  if {$index eq  2} {set ::highAff::CHAIN [split [lindex $argv $index] ,]}
  if {$index eq  3} {set ::highAff::RESID [lindex $argv $index]}
  if {$index eq  4} {set ::highAff::PROBE [split [lindex $argv $index] ,]} 
  if {$index eq 5} {set ::highAff::CUTOFF [lindex $argv $index]}
  if {$index eq 6} {set ::highAff::STEP [lindex $argv $index]}
}

if { $argc < 7 } {
  # 1st argument counted is the command itself so this should be one more than the number above
  puts "$index input arguments were provided on the command line."
  puts "The remaining values will be taken from the script."
}

if { $argc eq 2 && [lindex $argv 0] eq "help" } {
  puts ""
  puts "highaffresid.tcl takes up to 11 arguments, which must each contain no spaces. "
  puts "If you would like to include multiple values just put commas in between"
  puts "and the program will split them automatically."
  puts ""
  puts "The arguments are the following in a fixed order:"
  puts ""
  puts "struc: path to PDB file(s) for the starting structure of druggability simulation(s)"
  puts "you can provide either one path alone or one path per DCD file"
  puts "default value: ../md/sim.pdb"
  puts ""
  puts "dcd_in: path to DCD file(s) for druggability simulation(s)"
  puts "any number of DCD files can be provided as long as you meet the rule above about PDB files"
  puts "default value: ../md/sim.dcd"
  puts ""
  puts "CHAIN: the chains to analyze"
  puts "by default, all chains are analyzed"
  puts ""
  puts "RESID: the residue to analyze"
  puts "by default, all residues are analyzed"
  puts ""
  puts "PROBE: the probes to analyze"
  puts "by default, all probes are analyzed"
  puts ""
  puts "CUTOFF: the binding value cutoff for assigning high affinity residues"
  puts "default value is 500"
  puts ""
  exit
}

proc findHighAffResids { struc dcd_list CHAINS PROBES RESID } {

  set ofile   [open _out-$CHAIN-$PROBE.dat w]
  set ofile2  [open _out-detail-$CHAIN-$PROBE.dat w]
  set count 0
  set interval $STEP

  # start loop for dcd
  foreach dcd_in $dcd_list {

    mol load pdb $struc
    mol addfile $dcd_in type dcd first 0 last -1 step $STEP waitfor -1

    set argNit [atomselect top "resname $PROBE and not hydrogen"]
    set first_frame 0
    set num_frames [molinfo top get numframes] 

    set DISTSUM 0
    # start loop for frame
    for {set f $first_frame} {$f < $num_frames} {incr f} {
      incr count
      set time [expr ($count-1)*$interval]
      molinfo top set frame $f
      set gluOxy [atomselect top "chain $CHAIN and protein and resid $RESID and not hydrogen"]
      $gluOxy frame $f
      $argNit frame $f
      set MECO  [measure contacts $CUTOFF $gluOxy $argNit]
      set NNNN  [llength [lindex $MECO 0]] 
      if {$NNNN > 0} {  
        for {set mmm 0} {$mmm < $NNNN } {incr mmm} {
        set KKKA  [lindex [lindex $MECO 0 ] $mmm ]
        set KKKB  [lindex [lindex $MECO 1 ] $mmm ]
        set DIST  [measure bond [list [list $KKKA] [list $KKKB]]]    
        set DISTSUM [expr $DISTSUM+ 1/(($DIST)*($DIST))]

        set DISTV [expr 1/(($DIST)*($DIST))]
        set aato1 [expr $KKKA+1]
        set aato2 [expr $KKKB+1]
        set phe1 [atomselect [molinfo top get id] "serial $aato2"]
        set pha1 [atomselect [molinfo top get id] "serial $aato1"]
        set phe_vec1 [lindex [$phe1 get {resname}] 0]
        set phe_vec2 [lindex [$phe1 get {resid}] 0]
        set phe_vec3 [lindex [$phe1 get {name}] 0]
        set pha_vec1 [lindex [$pha1 get {resname}] 0]
        set pha_vec2 [lindex [$pha1 get {resid}] 0]
        set pha_vec3 [lindex [$pha1 get {name}] 0]
        puts $ofile2 "$f $aato1 $pha_vec1 $pha_vec2 $pha_vec3  $aato2 $phe_vec1 $phe_vec2 $phe_vec3  $DIST $DISTV"
        }
      }
    }
    # end loop for frame
    puts $ofile "$RESID  $DISTSUM"
    animate delete all
  }
  # end loop dcd list
  flush $ofile
  close $ofile
}
# end proc

if {[catch { 
  findHighAffResids $::highAff::struc $::highAff::dcd_in $::highAff::CHAIN $::highAff::PROBE $::highAff::RESID
} errvar ]} {
  puts $errvar
} else {
  puts "High affinity residue analysis completed successfully"
}
puts ""

exit
