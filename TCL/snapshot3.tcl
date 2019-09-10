# 2019-06-30 Bahar LAB, University of Pittsburgh  bahar@pitt.edu
# Written by Ji Young Lee
# Write PDB files for selected snapshots with dominant interactions from
# Druggability molecular dynamics simulations

set pdbfile APDB
set dcdfile ADCD

mol load pdb $pdbfile 
mol addfile  $dcdfile type dcd first 0 last -1 step SSTEP waitfor -1

animate write pdb pro-RRR-AAA.pdb beg AAA end AAA waitfor all sel [atomselect top "protein"]
animate write pdb lig-RRR-AAA-BBB-CCC.pdb beg AAA end AAA waitfor all sel [atomselect top "resname BBB and resid CCC"]

exit
