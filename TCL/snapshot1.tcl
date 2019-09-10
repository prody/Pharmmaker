# 2019-06-30 Bahar LAB, University of Pittsburgh  bahar@pitt.edu
# Written by Ji Young Lee
# Calculate interaction statistics from snapshots of 
# Druggability molecular dynamics simulations
# and rank dominant interactions

set ofile [open ligbo-RRR.dat w]
# Check cms
set struc      APDB
set dcd_list { ADCD}

#write a header for the datafile 
#  puts $ofile "#frame    Calpha_dist    N-O_dist    N-O_dist    N-O_dist    N-O_dist    N-O_dist    N-O_dist"

set count 0
# Check
set interval SSTEP 
foreach dcd_in $dcd_list {

# Check  step is related to interval
mol load pdb $struc
mol addfile $dcd_in type dcd first 0 last -1 step SSTEP waitfor -1

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

        set gluOxy [atomselect top "protein and chain CCC and resid BBB and not hydrogen"]
	$gluOxy frame $f
	set Olist [$gluOxy list]
############ Check distance cutoff
set argNit [atomselect top "resname AAA and not hydrogen"]
        $argNit frame $f
        set Nlist [$argNit list]
        set NNNN 0
	foreach atom1 $Olist {
		foreach atom2 $Nlist {
			set NOdist [measure bond [list [list $atom1] [list $atom2]]]
        if {$NOdist < CUTOFF} {
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
        animate write pdb v-RRR-$f.pdb beg $f end $f sel [atomselect top "protein or (resname AAA and resid $phe_vec2)"]
        puts $ofile "$f $aato1 $pha_vec1 $pha_vec2 $pha_vec3  $aato2 $phe_vec1 $phe_vec2 $phe_vec3  $NOdist"}
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
