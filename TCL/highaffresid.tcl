# 2019-06-30 Bahar LAB, University of Pittsburgh  bahar@pitt.edu
# Written by Ji Young Lee
# Calculate high affinity residues with trajectory of 
# Druggability molecular dynamics simulations

############ Check and adjust #################
# PDB file name and path
set struc      APDB                   
# DCD file name and path
set dcd_list { ADCD }
# protein chain ID in pdb
set CHAIN      ACHAIN
# residue number of protein for analysis
set RESID      ARESID
# Probe molecule ID in pdb
set PROBE      APROBE
# cutoff between residue and protein    
set CUTOFF     ACUTOFF
# frequency of dcd for analysis
set STEP       ASTEP
############ Check and adjust #################

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
exit
