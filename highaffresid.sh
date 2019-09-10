# 2019-06-30 Bahar LAB, University of Pittsburgh  bahar@pitt.edu
# Written by Ji Young Lee
# Calculate high affinity residues with trajectory of 
# Druggability molecular dynamics simulations

######### Check here  #######
# PDB file name and path
PDB='../drugui-simulation/protein-probe.pdb'
# DCD file name and path
DCD='../drugui-simulation/protein-probe.dcd'
# protein chain ID in pdb
CHAIN='A,B'    
# first residue number of protein for analysis
FIRSTRESID=5
# last residue number of protein for analysis
LASTRESID=261    # 261
# probe molecule ID in pdb
PROBE='IPRO,ACAM,ACTT,IPAM,IBTN,IMID' 
# cutoff between residue and probe      
CUTOFF=4.0
# frequency of dcd for analysis
STEP=1
######### Check here  #######

for FCHAIN in $CHAIN
do
  for FPROBE in $PROBE
  do
    for (( FF = $FIRSTRESID ; FF <= $LASTRESID  ; FF++ ))
    do
    env VMDARGS='text with blanks' vmd -dispdev text -e TCL/highaffresid.tcl $PDB $DCD $FCHAIN $FF $FPROBE $CUTOFF $STEP
    cat _out-$FCHAIN-$FPROBE.dat        >> out-$FCHAIN-$FPROBE.dat
    cat _out-detail-$FCHAIN-$FPROBE.dat >> out-detail-$FCHAIN-$FPROBE.dat
    done
  done
done

mkdir highaffresid
mv out-* highaffresid
