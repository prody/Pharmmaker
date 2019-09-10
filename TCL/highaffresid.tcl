# 2019-06-30 Bahar LAB, University of Pittsburgh  bahar@pitt.edu
# Written by Ji Young Lee
# Calculate high affinity residues with trajectory of 
# Druggability molecular dynamics simulations

namespace eval ::highAff:: {
  namespace export highAff

  ############ Check and adjust #################
  # PDB file name and path
  variable struc      '../drugui-simulation/protein-probe.pdb'
  # DCD file name and path
  variable dcd_list { '../drugui-simulation/protein-probe.dcd' }
  # protein chain ID in pdb
  variable CHAIN      all
  # first residue number of protein for analysis
  variable RESIDFIRST first
  # last residue number of protein for analysis
  variable RESIDLAST  last  
  # Probe molecule ID in pdb
  variable PROBE      all
  # cutoff between residue and protein    
  variable CUTOFF     4.0
  # frequency of dcd for analysis
  variable STEP       1
  ############ Check and adjust #################
}

# Take parameter values from input arguments as far as possible
for {set index 0} {$index < $argc -1} {incr index} {
  if {$index eq  0} {set ::highAff::struc [split [lindex $argv $index] ,]}
  if {$index eq  1} {set ::highAff::dcd_in [split [lindex $argv $index] ,]}
  if {$index eq  2} {set ::highAff::CHAIN [split [lindex $argv $index] ,]}
  if {$index eq  3} {set ::highAff::RESIDFIRST [split [lindex $argv $index] ,]}
  if {$index eq  4} {set ::highAff::RESIDLAST [split [lindex $argv $index] ,]}
  if {$index eq  5} {set ::highAff::PROBE [split [lindex $argv $index] ,]} 
  if {$index eq 6} {set ::highAff::CUTOFF [lindex $argv $index]}
  if {$index eq 7} {set ::highAff::STEP [lindex $argv $index]}
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
  puts "RESIDFIRST: the first residue to include in the analysis for each chain"
  puts "by default, the first residue in each chain is taken"
  puts "If you provide just one value, it will be used for all chains"
  puts ""
  puts "RESIDLAST: the last residue to include in the analysis for each chain"
  puts "by default, the last residue in each chain is taken"
  puts "If you provide just one value, it will be used for all chains"
  puts ""
  puts "PROBE: the probes to analyze"
  puts "by default, all probes are analyzed"
  puts ""
  puts "CUTOFF: the binding value cutoff for assigning high affinity residues"
  puts "default value is 500"
  puts ""
  exit
}

proc findHighAffResids { struc dcd_list CHAINS PROBES RESID STEP } {

  # start loop for chains
  set chainNum 0
  foreach CHAIN $CHAINS {
    puts ""
    puts "starting loop for chain $CHAIN"
    puts ""

    # start loop for probes
    foreach PROBE $PROBES {
      set probeSel [atomselect top "resname $PROBE and not hydrogen"]
      $probeSel frame $f

      puts "This analysis will run from residue [lindex $RESIDFIRST $chainNum] to residue [lindex $RESIDLAST $chainNum]"
      puts ""
      set resids {}
      for {set rrr [lindex $RESIDFIRST $chainNum] } {$rrr <= [lindex $RESIDLAST $chainNum] } {incr rrr} {
        lappend resids $rrr
      }

      # start loop for residue
      foreach rrr $resids {
        set DISTSUM 0
        set ofile   [open _out-$CHAIN-$PROBE.dat w]
        set ofile2  [open _out-detail-$CHAIN-$PROBE.dat w]
        set count 0

        # start loop for dcd
        foreach dcd_in $dcd_list {

          mol load pdb $struc
          mol addfile $dcd_in type dcd first 0 last -1 step $STEP waitfor -1

          set first_frame 0
          set num_frames [molinfo top get numframes] 

          # start loop for frame
          for {set f $first_frame} {$f < $num_frames} {incr f} {
            incr count
            set time [expr ($count-1)*$STEP]
            molinfo top set frame $f

            set residSel [atomselect top "chain $CHAIN and protein and resid $RESID and not hydrogen"]
            $residSel frame $f

            set MECO  [measure contacts $CUTOFF $residSel $probeSel]
            set numConts  [llength [lindex $MECO 0]] 
            if {$numConts > 0} {  
              for {set mmm 0} {$mmm < $numConts } {incr mmm} {
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
        }
        # end loop dcd list

        $residSel delete

        puts $ofile "$RESID  $DISTSUM"
        animate delete all
        flush $ofile
        close $ofile
      }
      # end loop for residue
      $probeSel delete

      if { $chainNum eq 0 } {
        set pfile [open probe-list.dat a]
        puts $pfile $PROBE
        flush $pfile
        close $pfile
      }
    }
    # end loop for probes
    incr chainNum

    set cfile [open chain-list.dat a]
    puts $cfile $CHAIN
    flush $cfile
    close $cfile
  }
  # end loop for chains

}
# end proc

##############################################################################
# Declare some more procs for automatically finding chains, probes and residues 

proc findChainIDs { struc } {
  mol load pdb $struc
  incr ::highAff::mol_ID
  set sel [atomselect top all]
  set chains [lsort -ascii -unique [$sel get chain]] 
  set ::highAff::CHAIN [lrange $chains 0 1]
  $sel delete
  animate delete all
  mol delete top
}

proc findProbeTypes { struc } {
  mol load pdb $struc
  incr ::highAff::mol_ID
  set sel [atomselect top "chain P"]
  set ::highAff::PROBE [lsort -ascii -unique [$sel get resname]] 
  $sel delete
  animate delete all
  mol delete top
}

proc findResidRange { struc chains } {
  mol load pdb $struc
  incr ::highAff::mol_ID
  
  # Loop through the chains and append to the list if needed
  # To know if we need it, we append without overwriting first

  foreach chain $chains {
    set sel [atomselect top "chain $chain and name CA"]
    set resids [$sel get resid]

    if { [lindex $::highAff::RESIDFIRST 0] eq "first" } {
      lappend ::highAff::RESIDFIRST [lindex $resids 0]
    }

    if { [lindex $::highAff::RESIDLAST 0] eq "last" } {
      lappend ::highAff::RESIDLAST [lindex $resids end]
    }
  }

  # Remove the 0th element if necessary

  if { [lindex $::highAff::RESIDFIRST 0] eq "first" } {
    set ::highAff::RESIDFIRST [lrange $::highAff::RESIDFIRST 1 end]
  }

  if { [lindex $::highAff::RESIDLAST 0] eq "last" } {
    set ::highAff::RESIDLAST [lrange $::highAff::RESIDLAST 1 end]
  }

  $sel delete
  animate delete all
  mol delete top
}

##############################################################################
# Call those procs if values aren't already provided

if { $::highAff::CHAIN eq "all" } {
  findChainIDs $::highAff::struc
}

if { $::highAff::PROBE eq "all" } {
  findProbeTypes $::highAff::struc
}

if { $::highAff::RESIDFIRST eq "first" || $::highAff::RESIDLAST eq "last" } {
  findResidRange $::highAff::struc $::highAff::CHAIN
}

##############################################################################
# Call the the main proc with error handling

if {[catch { 
  findHighAffResids $::highAff::struc $::highAff::dcd_in $::highAff::CHAIN $::highAff::PROBE $::highAff::RESID $::highAff::STEP
} errvar ]} {
  puts $errvar
} else {
  puts "High affinity residue analysis completed successfully"
}
puts ""

exit
