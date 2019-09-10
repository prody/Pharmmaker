# 2019-06-30 Bahar LAB, University of Pittsburgh  bahar@pitt.edu
# Written by Ji Young Lee
# Calculate high affinity residues with trajectory of 
# Druggability molecular dynamics simulations

namespace eval ::highAff:: {
  namespace export highAff

  ############ Check and adjust #################
  #  0. PDB file name and path
  variable struc      '../drugui-simulation/protein-probe.pdb'
  #  1. DCD file name and path
  variable dcd_list { '../drugui-simulation/protein-probe.dcd' }
  #  2. binding value cutoff for assigning high affinity residues
  variable bv_cutoff 500
  #  3. protein chain ID in pdb
  variable CHAIN      all
  #  4. first residue number of protein for analysis
  variable RESIDFIRST first
  #  5. last residue number of protein for analysis
  variable RESIDLAST  last  
  #  6. Probe molecule ID in pdb
  variable PROBE      all
  #  7. cutoff between residue and protein    
  variable CUTOFF     4.0
  #  8. frequency of dcd for analysis
  variable STEP       1
  #  9. first frame number for analysis
  variable frameFIRST first
  # 10. last frame number for analysis
  variable frameLAST  last 
  ############ Check and adjust #################
}

# Take parameter values from input arguments as far as possible
for {set index 0} {$index < $argc -1} {incr index} {
  if {$index eq  0} {set ::highAff::struc [split [lindex $argv $index] ,]}
  if {$index eq  1} {set ::highAff::dcd_in [split [lindex $argv $index] ,]}
  if {$index eq  2} {set ::highAff::bv_cutoff [split [lindex $argv $index] ,]}
  if {$index eq  3} {set ::highAff::CHAIN [split [lindex $argv $index] ,]}
  if {$index eq  4} {set ::highAff::RESIDFIRST [split [lindex $argv $index] ,]}
  if {$index eq  5} {set ::highAff::RESIDLAST [split [lindex $argv $index] ,]}
  if {$index eq  6} {set ::highAff::PROBE [split [lindex $argv $index] ,]} 
  if {$index eq  7} {set ::highAff::CUTOFF [lindex $argv $index]}
  if {$index eq  8} {set ::highAff::STEP [lindex $argv $index]}
  if {$index eq  9} {set ::highAff::frameFIRST [split [lindex $argv $index] ,]}
  if {$index eq 10} {set ::highAff::frameLAST [split [lindex $argv $index] ,]}  
}

if { $argc < 12 } {
  # 1st argument is the command then we count from 0 so this should be two more than the number above
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
  puts "bv_cutoff: the binding value cutoff for assigning high affinity residues"
  puts "default value is 500"
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
  puts "frameFIRST: the first frame to include in the analysis"
  puts "by default, the first frame is taken"
  puts ""
  puts "frameLAST: the last frame to include in the analysis"
  puts "by default, the last frame is taken"
  puts ""
  exit
}

proc findHighAffResids { struc dcd_list CHAINS PROBES RESIDFIRST RESIDLAST STEP BV_CUTOFF CUTOFF frameFIRST frameLAST } {

  if {[file exists highaffresids] eq 0} {
    file mkdir highaffresids
  }

  set dcdNum 0
  # start loop for dcd
  foreach dcd_in $dcd_list {

    mol load pdb $struc
    mol addfile $dcd_in type dcd first 0 last -1 step $STEP waitfor -1

    # start loop for frame
    for {set f $frameFIRST } {$f <= $frameLAST } {set f [expr {$f + $STEP}]} {
      molinfo top set frame $f
      if { [expr $f % 100] eq 0 } {
        puts "starting loop for frame $f"
      }

      # start loop for chains
      set chainNum 0
      foreach CHAIN $CHAINS {
        #puts ""
        #puts "starting loop for chain $CHAIN"
        #puts ""

        # start loop for probes
        set probeNum 0
        foreach PROBE $PROBES {
          set probeSel [atomselect top "resname $PROBE and not hydrogen"]
          $probeSel frame $f

          #puts "This analysis will run from residue [lindex $RESIDFIRST $chainNum] to residue [lindex $RESIDLAST $chainNum]"
          #puts ""

          set highaffresids {}
          # start loop for residue
          set resCounter 0
          for {set rrr [lindex $RESIDFIRST $chainNum] } {$rrr <= [lindex $RESIDLAST $chainNum] } {incr rrr} {

            set residSel [atomselect top "chain $CHAIN and protein and resid $rrr and not hydrogen"]
            $residSel frame $f

            if { $f eq $frameFIRST && $dcdNum eq 0 } {
              # We use a multi-dimensional associative array for DISTSUM values.
              # We index it by counting chains, probes and residues
              # These are counted with the variables chainNum, probeNum and resCounter
              # Each DISTSUM value should be the sum of inverse square distances within a cutoff (binding values) 
              set DISTSUMS($chainNum,$probeNum,$resCounter) 0
              # It gets initialised the first time we fill it
              # We then update DISTSUM from it each time and fill it back up again
            }
            set DISTSUM $DISTSUMS($chainNum,$probeNum,$resCounter)

            set MECO  [measure contacts $CUTOFF $residSel $probeSel]
            set numConts  [llength [lindex $MECO 0]] 
            if {$numConts > 0} {  
              for {set mmm 0} {$mmm < $numConts } {incr mmm} {
              set KKKA  [lindex [lindex $MECO 0 ] $mmm ]
              set KKKB  [lindex [lindex $MECO 1 ] $mmm ]
              set DIST  [measure bond [list [list $KKKA] [list $KKKB]]]    
              set DISTSUM [expr $DISTSUM + 1/(($DIST)*($DIST))]

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
              }
            }

            set DISTSUMS($chainNum,$probeNum,$resCounter) $DISTSUM

            $residSel delete

            if { $f eq $frameLAST } {
              set ofile   [open highaffresids/$CHAIN-$PROBE.dat a]
              puts $ofile "$rrr  $DISTSUM"
              flush $ofile
              close $ofile

              if { $DISTSUM > $BV_CUTOFF } {
                lappend highaffresids $rrr
              }
            }

            incr resCounter
          }
          # end loop for residue
          $probeSel delete

          if { $chainNum eq 0 && $f eq $frameFIRST && $dcdNum eq 0 } {
            set pfile [open probe-list.dat a]
            puts $pfile $PROBE
            flush $pfile
            close $pfile
          }

          if {[llength $highaffresids] > 0} {
            set hfile  [open highaffresids/$CHAIN-$PROBE-highaffresid.dat a]
            puts $hfile "$highaffresids"
            flush $hfile
            close $hfile

            set hfile2  [open highaffresid.dat a]
            puts $hfile2 "$PROBE $CHAIN $highaffresids"
            flush $hfile2
            close $hfile2
          }

          incr probeNum
        }
        # end loop for probes
        incr chainNum

        if { $f eq $frameFIRST && $dcdNum eq 0 } {
          set cfile [open chain-list.dat a]
          puts $cfile $CHAIN
          flush $cfile
          close $cfile
        }
      }
      # end loop for chains

    }
    # end loop for frame
  }
  # end loop dcd list
  animate delete all
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

proc findFrameRange { struc dcd_in } {
  mol load pdb $struc
  incr ::highAff::mol_ID
  mol addfile $dcd_in type dcd first 0 last -1 step 1 waitfor -1

  if { $::highAff::frameFIRST eq "first" } {
    set ::highAff::frameFIRST 1
  }
    
  if { $::highAff::frameLAST eq "last" } {
    set num_frames [molinfo top get numframes]
    set ::highAff::frameLAST [expr $num_frames - 1]
  }

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

if { [llength $::highAff::RESIDFIRST] eq 1 } {
  set RESIDFIRST $::highAff::RESIDFIRST
  set ::highAff::RESIDFIRST {}
  foreach CHAIN $::highAff::CHAIN {
    lappend ::highAff::RESIDFIRST $RESIDFIRST
  }
}

if { [llength $::highAff::RESIDLAST] eq 1 } {
  set RESIDLAST $::highAff::RESIDLAST
  set ::highAff::RESIDLAST {}
  foreach CHAIN $::highAff::CHAIN {
    lappend ::highAff::RESIDLAST $RESIDLAST
  }
}

if { $::highAff::frameFIRST eq "first" || $::highAff::frameLAST eq "last" } {
  findFrameRange $::highAff::struc $::highAff::dcd_in
}

##############################################################################
# Call the the main proc with error handling

if {[catch { 
  findHighAffResids $::highAff::struc $::highAff::dcd_in $::highAff::CHAIN $::highAff::PROBE $::highAff::RESIDFIRST $::highAff::RESIDLAST $::highAff::STEP $::highAff::bv_cutoff $::highAff::CUTOFF $::highAff::frameFIRST $::highAff::frameLAST
} errvar ]} {
  puts $errvar
} else {
  puts "High affinity residue analysis completed successfully"
}
puts ""

exit
