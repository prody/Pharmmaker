# 2019-06-30 Bahar LAB, University of Pittsburgh  bahar@pitt.edu
# Written by Ji Young Lee
# Select snapshots with many of the most dominant interactions from
# Druggability molecular dynamics simulations

set ofile [open _ligbo-ok.dat w]
# Check cms
set input    ./v-com-ok.pdb

set count 0
set interval SSTEP 

animate read pdb $input beg 0 end -1 skip SSTEP waitfor all

set first_frame 0
set num_frames [molinfo top get numframes] 

for {set f $first_frame} {$f < $num_frames} {incr f} {
	incr count
	set time [expr ($count-1)*$interval]
	molinfo top set frame $f

        set gluOxy [atomselect top "AAA"]
	$gluOxy frame $f
	set Olist [$gluOxy list]
############ Check distance cutoff
set argNit [atomselect top "BBB"]
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
#        animate write pdb xok-$f.pdb beg $f end $f sel [atomselect top "all"]
        puts $ofile "$f $aato1 $pha_vec1 $pha_vec2 $pha_vec3  $aato2 $phe_vec1 $phe_vec2 $phe_vec3  $NOdist"}
		}
	}	
############
}
#end loop over frames
animate delete all
#end dtr list
flush $ofile
close $ofile
exit
