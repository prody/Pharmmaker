# 2019-06-30 Bahar LAB, University of Pittsburgh  bahar@pitt.edu
# Written by Ji Young Lee, and edited by James Michael Krieger
# Calculate high affinity residues from trajectories of 
# Druggability molecular dynamics simulations

package require namdenergy
package require readcharmmpar
package require readcharmmtop

namespace eval ::highAff:: {
  namespace export highAff

  ############ Default parameter values (change or provide) #################
  # 0. PDB file name and path
  variable struc  ../md/sim.pdb     
  # 1. DCD file name and path
  variable dcd_in ../md/sim.dcd
  # 2. resids or a directory containing PDB files for refining a region of interest
  variable region_refiner None
  # 3. protein chain ID in pdb
  variable CHAIN all
  # 4. Probe molecule ID in pdb
  variable PROBE all
  # 5. binding value cutoff for assigning high affinity residues
  variable bv_cutoff 500
  # 6. cutoff for selecting residues as being near hotspots
  variable region_point_residue_cutoff 8
  # 7. first residue number of protein for analysis
  variable RESIDFIRST first
  # 8. last residue number of protein for analysis
  variable RESIDLAST  last  
  # 9. first frame number for analysis
  variable frameFIRST first
  # 10. last frame number for analysis
  variable frameLAST  last  
  ############ Default parameter values (change or provide) #################

  # namespace variable for storing residue range generated from RESIDFIRST, RESIDLAST and/or hotspots
  variable resids

  # namespace variable for storing frame range generated from frameFIRST and frameLAST
  variable frameRange

  # mol ID for top molecule, which we need to track globally for making selections from non-top molecules
  variable mol_ID
}


# Take parameter values from input arguments as far as possible
for {set index 0} {$index < $argc -1} {incr index} {
  if {$index eq  0} {set ::highAff::struc [split [lindex $argv $index] ,]}
  if {$index eq  1} {set ::highAff::dcd_in [split [lindex $argv $index] ,]}
  if {$index eq  2} {set ::highAff::region_refiner [split [lindex $argv $index] ,]}
  if {$index eq  3} {set ::highAff::CHAIN [split [lindex $argv $index] ,]}
  if {$index eq  4} {set ::highAff::PROBE [split [lindex $argv $index] ,]}
  if {$index eq  5} {set ::highAff::bv_cutoff [split [lindex $argv $index] ,]}
  if {$index eq  6} {set ::highAff::region_point_residue_cutoff [lindex $argv $index]}
  if {$index eq  7} {set ::highAff::RESIDFIRST [split [lindex $argv $index] ,]}
  if {$index eq  8} {set ::highAff::RESIDLAST [split [lindex $argv $index] ,]}
  if {$index eq  9} {set ::highAff::frameFIRST [split [lindex $argv $index] ,]}
  if {$index eq 10} {set ::highAff::frameLAST [split [lindex $argv $index] ,]}  
}

if {$index eq 8} {set ::highAff::RESIDLAST $::highAff::RESIDFIRST}

set refiners {}
foreach refiner $::highAff::region_refiner {
  if {[catch { 
    set pdbfiles [glob -directory $refiner *site_?_soln_?.pdb]
    foreach pdbfile $pdbfiles {
      lappend refiners $pdbfile
    }
  } errvar ]} {
    lappend refiners $refiner
  }
}
set ::highAff::region_refiner $refiners

if { $argc < 12 } {
  # 1st argument counted is the command itself
  puts "$index input arguments were provided."
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
  puts "region_refiner: resids or a directory containing PDB files for refining a region of interest"
  puts "default value is None"
  puts ""
  puts "CHAIN: the chains to analyze"
  puts "by default, all chains are analyzed"
  puts ""
  puts "PROBE: the probes to analyze"
  puts "by default, all probes are analyzed"
  puts ""
  puts "bv_cutoff: the binding value cutoff for assigning high affinity residues"
  puts "default value is 500"
  puts ""
  puts "region_point_residue_cutoff: the cutoff for interactions of region-defining points with residues in the initial structure"
  puts "for region of interest selection from hotspots. Default value is 8 A."
  puts ""
  puts "RESIDFIRST: the first residue to include in the analysis for each chain"
  puts "by default, the first residue in each chain is taken"
  puts "If you provide just one value, it will be used for all chains"
  puts ""
  puts "RESIDLAST: the last residue to include in the analysis for each chain"
  puts "by default, the last residue in each chain is taken"
  puts "If you provide just one value, it will be used for all chains"
  puts ""
  puts "To set a region of interest using resids, provide the same list to RESIDFIRST and RESIDLAST or leave RESIDLAST blank"
  puts ""
  puts "frameFIRST: the first frame to include in the analysis"
  puts "by default, the first frame is taken"
  puts ""
  puts "frameLAST: the last frame to include in the analysis"
  puts "by default, the last frame is taken"
  puts ""
  puts "To set a frame range of interest, provide the same list to frameFIRST and frameLAST"
  puts ""
  exit
}

proc findHighAffResids { sites strucs dcd_list CHAINS PROBES RESIDFIRST RESIDLAST cutoff frameFIRST frameLAST } {

  set bestSite 0
  set bestScore 0
  set siteNum 0
  # start loop for sites
  foreach site $sites {
    set siteScore 0

    puts ""
    puts "starting loop for site $site"

    # start loop for chains
    set chainNum 0
    foreach CHAIN $CHAINS {
      puts ""
      puts "starting loop for chain $CHAIN"
      puts ""

      # make a list of resids for the resids foreach loop and handle regions of interest
      if { $RESIDFIRST eq $RESIDLAST } {
        puts "A region of interest was selected by providing a list of residues. The same list will be used for all chains."
        puts "The residues of interest are $RESIDFIRST (a total of [llength $RESIDFIRST] residues)"
        # This only works with running the tcl script on one chain at a time.
        set ::highAff::resids $RESIDFIRST
      } else {
        puts "This analysis will run from residue [lindex $RESIDFIRST $chainNum] to residue [lindex $RESIDLAST $chainNum]"
        puts ""
        set ::highAff::resids {}
        for {set rrr [lindex $RESIDFIRST $chainNum] } {$rrr <= [lindex $RESIDLAST $chainNum] } {incr rrr} {
          lappend ::highAff::resids $rrr
        }
      }

      # make a list of frames for the frames foreach loop
      if { $frameFIRST eq $frameLAST } {
        puts "A frame range was selected by providing a list of frames."
        puts "The frames of interest are $frameFIRST (a total of [llength $frameFIRST] frames)"
        # This only works with running the tcl script on one chain at a time.
        set ::highAff::frameRange $frameFIRST
      } else {
        puts "This analysis will run from frame $frameFIRST to frame $frameLAST"
        puts ""
        set ::highAff::frameRange {}
        for {set rrr $frameFIRST } {$rrr <= $frameLAST } {incr rrr} {
          lappend ::highAff::frameRange $rrr
        }
      }

      refineResidRange $::highAff::resids $CHAIN $site
      puts ""
      puts "The residues selected as the region of interest in chain $CHAIN for site $site are $::highAff::resids (a total of [llength $::highAff::resids] residues)"
      puts ""

      set siteList1 [split $site "/"]
      set siteString1 [lindex $siteList1 end]
      set siteList2 [split $siteString1 "."]
      set site_string [lindex $siteList2 0]

      if {[string tolower $site_string]=="none"} {
        set site_string "highaffresid"
      }

      if {[file exists $site_string] eq 0} {
        file mkdir $site_string
        exec ln -s ../$site $site_string/
      }
      set regionFile [open "$site_string/region-residues.dat" a]
      puts $regionFile "$CHAIN $::highAff::resids"
      flush $regionFile
      close $regionFile

      set regionFile2 [open "region-residues.dat" a]
      puts $regionFile2 "$site_string $CHAIN $::highAff::resids"
      flush $regionFile2
      close $regionFile2

      if { [llength $::highAff::resids] eq 0 } {
        incr chainNum
        # continue on to next chain
        continue
      }

      # start loop for dcd
      set dcdNum 0
      foreach dcd_in $dcd_list {

        if {$siteNum eq 0} {
          set dfile [open "dcd-list.dat" a]
          puts $dfile $dcd_in
          flush $dfile
          close $dfile

          if { [llength $strucs] > 1 } {
            set struc [lindex $strucs $dcdNum]
          } else {
            set struc $strucs
          }

          if {[expr { $dcdNum < [llength $strucs] }] && $chainNum eq 0} {
            set sfile [open "struc-list.dat" a]
            puts $sfile $struc
            flush $sfile
            close $sfile
          }
        }

        mol load pdb $struc
        incr ::highAff::mol_ID
        mol addfile $dcd_in type dcd first 0 last -1 step 1 waitfor -1

        set first_frame [lindex $::highAff::frameRange 0]
        set num_frames [llength $::highAff::frameRange]

        puts ""
        # start loop for frame
        foreach f $::highAff::frameRange {
          molinfo top set frame $f
          #puts ""
          if { [expr $f % 100] eq 0} {
            puts "starting loop for frame $f of $num_frames (using chain $CHAIN residues near site $site)"
          }

          # start loop for probes
          set probeNum 0
          foreach PROBE $PROBES {
            #puts ""
            #puts "starting loop for probe $PROBE"
            set highaffresids {}

            set argNit [atomselect top "resname $PROBE and not hydrogen"]
            $argNit frame $f
            
            set resCounter 0
            # start loop for residue
            foreach rrr $::highAff::resids {
              #puts "starting loop for residue $rrr in chain $CHAIN with probe $PROBE"

              if { $f eq $first_frame } {
                # We use a multi-dimensional associative array for DISTSUM values.
                # We index it by counting sites, chains, probes and residues
                # These are counted with the variables siteNum, chainNum, probeNum and resCounter
                # Each DISTSUM value should be the sum of inverse square distances within a cutoff (binding values) 
                set DISTSUMS($siteNum,$chainNum,$probeNum,$resCounter) 0
                # It gets initialised the first time we fill it
                # We will go through this multi-dimensional list to write files later
                # We will also decide upon the best site for further analysis at that stage

                # We likewise create a multi-dimensional associative array for residue identifiers
                # We index it as above
                set RESIDS($siteNum,$chainNum,$probeNum,$resCounter) $rrr
              }
              set DISTSUM $DISTSUMS($siteNum,$chainNum,$probeNum,$resCounter)

              set gluOxy [atomselect top "chain $CHAIN and protein and resid $rrr and not hydrogen"]
              $gluOxy frame $f

              set MECO  [measure contacts 4 $gluOxy $argNit]
              #puts $MECO
              set NNNN  [llength [lindex $MECO 0]] 
              if {$NNNN > 0} {  
                for {set mmm 0} {$mmm < $NNNN } {incr mmm} {
                  set KKKA  [lindex [lindex $MECO 0 ] $mmm ]
                  set KKKB  [lindex [lindex $MECO 1 ] $mmm ]
                  set DIST  [measure bond [list [list $KKKA] [list $KKKB]]]    
                  set DISTSUM [expr $DISTSUM + 1/(($DIST)*($DIST))]
                }
              }
              $gluOxy delete

              #puts $DISTSUM
              set DISTSUMS($siteNum,$chainNum,$probeNum,$resCounter) $DISTSUM
              if { $f eq [expr $num_frames - 1]} {
                #puts "$siteNum $site $chainNum $CHAIN $probeNum $PROBE $resCounter $rrr $DISTSUM"

                set ofile [open $site_string/$CHAIN-$PROBE.dat a]
                
                puts $ofile "$rrr  [format {%0.2f} $DISTSUM]"
                if { $DISTSUM > $cutoff } {
                  lappend highaffresids $rrr
                  set siteScore [ expr $siteScore + $DISTSUM ]
                }
                flush $ofile
                close $ofile

              }
              incr resCounter
            }
            # end loop for residue
            $argNit delete
            incr probeNum

            if { $f eq $first_frame && $siteNum eq 0 && $chainNum eq 0 } {
              set pfile [open probe-list.dat a]
              puts $pfile $PROBE
              flush $pfile
              close $pfile
            }

            if { $f eq [expr $num_frames - 1]} {

              set siteList1 [split $site "/"]
              set siteString1 [lindex $siteList1 end]
              set siteList2 [split $siteString1 "."]
              set site_string [lindex $siteList2 0]

              if {[llength $highaffresids] > 0} {
                set ofile2  [open highaffresid.dat a]
                puts $ofile2 "$site $PROBE $CHAIN $highaffresids"
                flush $ofile2
                close $ofile2
              }

              set ofile2  [open $site_string/$CHAIN-$PROBE-highaffresid.dat a]
              puts $ofile2 "$highaffresids"
              flush $ofile2
              close $ofile2

              if {[llength $highaffresids] > 0} {
                set ofile2  [open $site_string/highaffresid.dat a]
                puts $ofile2 "$PROBE $CHAIN $highaffresids"
                flush $ofile2
                close $ofile2
              }
            }
          }
          # end loop for probes

        }
        # end loop for frames

        animate delete all
        incr dcdNum
      }
      # end loop for dcds

      incr chainNum

      if { $siteNum eq 0 } {
        set cfile [open chain-list.dat a]
        puts $cfile $CHAIN
        flush $cfile
        close $cfile
      }
    }
    # end loop for chains

    set sfile [open site-scores.dat a]
    puts $sfile "[format {%0.2f} $siteScore] $site"
    flush $sfile
    close $sfile

    if { $siteScore > $bestScore } {
      set bestScore $siteScore
      set bestSite $site
    }

    incr siteNum
  }
  # end loop for sites

  puts ""
  puts "The best site is $bestSite"
  puts ""

  set siteList1 [split $bestSite "/"]
  set siteString1 [lindex $siteList1 end]
  set siteList2 [split $siteString1 "."]
  set site_string [lindex $siteList2 0]

  # Make symbolic links for directory and original PDB for best site
  exec ln -s $site_string best-site
}
# end proc


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

proc refineResidRange { resid_range chain refiner } {
  if {[string tolower $refiner]=="none"} {
    puts "Region refinment is not carried out"
  } else {

    mol load pdb $::highAff::struc
    incr ::highAff::mol_ID
    set pro [atomselect top "chain $chain"]

    foreach region_refiner $refiner {
    
      if {[catch { 
        mol load pdb $region_refiner
        #incr ::highAff::mol_ID later to access protein mol_ID more easily
        set rp_set [atomselect top "all"]
        set rp_resids [$rp_set get resid]
        $rp_set delete
        set rp_definer "file"
      } errvar ]} {
        puts "The region point refiner was a string not file"
        set rp_definer "string"
        set rp_resids $region_refiner
      }

      puts ""
      puts "region-defining point resids are $rp_resids"
  
      set ::highAff::resids {}
      foreach resid $resid_range {
        set residue [atomselect $::highAff::mol_ID "chain $chain and resid $resid"]

        foreach rp_resid $rp_resids {
          set rp [atomselect top "resid $rp_resid"]
          set conts [measure contacts $::highAff::region_point_residue_cutoff $rp $residue]
          set num_conts [llength [lindex $conts 0]]
          if {$num_conts > 0 && [lsearch $::highAff::resids $resid] < 0} {
            lappend ::highAff::resids $resid
          }
          $rp delete
        }
        $residue delete
      }

      if {$rp_definer=="file"} {
        mol delete top 
      }

      incr ::highAff::mol_ID
    }
    $pro delete
    mol delete top
  }
}

set ::highAff::mol_ID -1
# This means that we get to mol ID 0 when we increment it the first time

if { $::highAff::CHAIN eq "all" } {
  findChainIDs $::highAff::struc
}

if { $::highAff::PROBE eq "all" } {
  findProbeTypes $::highAff::struc
}

if { $::highAff::RESIDFIRST eq "first" || $::highAff::RESIDLAST eq "last" } {
  findResidRange $::highAff::struc $::highAff::CHAIN
}

if { $::highAff::frameFIRST eq "first" || $::highAff::frameLAST eq "last" } {
  findFrameRange $::highAff::struc $::highAff::dcd_in
}

set numChains [llength $::highAff::CHAIN]

if { [llength $::highAff::RESIDFIRST] eq 1 && $numChains > 1} {
  for {set index 1} {$index < $numChains} {incr index} {
    lappend ::highAff::RESIDFIRST [lindex $::highAff::RESIDFIRST 0]
  }
}

if { [llength $::highAff::RESIDLAST] eq 1 && $numChains > 1} {
  for {set index 1} {$index < $numChains} {incr index} {
    lappend ::highAff::RESIDLAST [lindex $::highAff::RESIDLAST 0]
  }
}

#puts $::highAff::struc
#puts $::highAff::dcd_in
#puts $::highAff::CHAIN
#puts $::highAff::PROBE
#puts $::highAff::RESIDFIRST
#puts $::highAff::RESIDLAST

# set variables for testing
######################################
#set strucs $::highAff::struc
#set dcd_list $::highAff::dcd_in
#set CHAINS $::highAff::CHAIN
#set PROBES $::highAff::PROBE
#set RESIDFIRST $::highAff::RESIDFIRST
#set RESIDLAST $::highAff::RESIDLAST
#set cutoff $::highAff::bv_cutoff
#set chainNum 0
#set CHAIN [lindex $CHAINS 0]
######################################

if {[catch { 
  findHighAffResids $::highAff::region_refiner $::highAff::struc $::highAff::dcd_in $::highAff::CHAIN $::highAff::PROBE $::highAff::RESIDFIRST $::highAff::RESIDLAST $::highAff::bv_cutoff $::highAff::frameFIRST $::highAff::frameLAST
} errvar ]} {
  puts $errvar
} else {
  puts "High affinity residue analysis completed successfully"
}
puts ""
exit
