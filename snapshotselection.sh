# 2019-06-30 Bahar LAB, University of Pittsburgh  bahar@pitt.edu
# Written by Ji Young Lee
# Select snapshots with many of the most dominant interactions from
# Druggability molecular dynamics simulations

######### Check here  #######
# PDB file name and path
PDB='\/data\/z1906-tool2020\/tutorial-new\/drugui-simulationRRR\/protein-probe.pdb'
# DCD file name and path RRR is variable
DCD='\/data\/z1906-tool2020\/tutorial-new\/drugui-simulationRRR\/protein-probe.dcd'
# DCD trajectory number for RRR
TRJNUM=' 1 '
# frame number for each dcd         
FRANUM=10000
# protein chain ID in pdb
#CHAIN=' A B '
# probe molecule ID in pdb
#PROBE=' IPRO ACAM ACTT IPAM IBTN IMID '
# cutoff between residue and probe      
#CUTOFF=4.0
# cutoff between probe and hot spot
#CUTOFF2=1.5
# frequency of dcd for analysis
STEP=1
# output directory             
OUTDIR=' s1 '
# cutoff of score for collecting snapshots
CUTOFF3=1000
######### Check here  #######

for FOUTDIR in $OUTDIR
do

#########
cat > _molexe  <<EOF
awk '{if (\$1 >= $CUTOFF3) print }' snapshot-$FOUTDIR/zlist-count  > ___test 
EOF
chmod 777 _molexe
./_molexe
rm _molexe

mkdir CUT$CUTOFF3

TNUM=`wc -l ___test | awk '{print $1}'`
for (( mm = 1 ; mm <= $TNUM  ; mm++ ))
do
sed -n -e "$mm,$mm p" ___test > ___test2
sed -e "s/outfr.dat/  /g" ___test2 > ___test3

FF=`cat ___test3 | awk '{print $2}'`
cp -r $FF CUT$CUTOFF3
done

rm ___t*
############

######### zlist-frame-c1000 zlist-frame-c1000-detail
for KK in $TRJNUM
do
for (( nn = 0   ; nn <= $FRANUM   ; nn++ ))
do

grep "SIM# $KK FRAME# $nn "  CUT$CUTOFF3/z.*/outfr.dat > ___test

CRAT=`wc -l ___test | awk '{print $1}'`

if [ $CRAT -eq $TNUM ]; then
cat   ___test >> ___test2
echo "SIM# $KK FRAME# $nn   $CRAT"  >> ___test3
fi

done
done

cp ___test2    snapshot-$FOUTDIR/zlist-frame-c$CUTOFF3-detail
cp ___test3    snapshot-$FOUTDIR/zlist-frame-c$CUTOFF3

rm  __*
\rm -r CUT$CUTOFF3
######### zlist-frame-c1000 zlist-frame-c1000-detail


# Extract protein and ligand ##

SNAPSNUM=`wc -l snapshot-$FOUTDIR/zlist-frame-c$CUTOFF3 | awk '{print $1}'`

for (( yy = 1 ; yy <= $SNAPSNUM  ; yy++ ))
do

sed -n -e "$yy,$yy p" snapshot-$FOUTDIR/zlist-frame-c$CUTOFF3  > _____tt
FF=`cat _____tt | awk '{print $2}'`
FFRAM=`cat _____tt | awk '{print $4}'`
grep "$FF FRAME# $FFRAM" snapshot-$FOUTDIR/zlist-frame-c$CUTOFF3-detail > _____tt2

TTTTNUM=`wc -l _____tt2 | awk '{print $1}'`

if [ $TTTTNUM -ge 1 ]; then
for (( bb = 1 ; bb <= $TTTTNUM  ; bb++ ))
do
sed -n -e  "$bb,$bb p" _____tt2 > _____tt3

PPROB=`cat _____tt3 | awk '{print $6}'`
PPRON=`cat _____tt3 | awk '{print $8}'`
sed -e "s/APDB/$PDB/g" -e "s/ADCD/$DCD/g" -e "s/SSTEP/$STEP/g" -e "s/AAA/$FFRAM/g" -e "s/BBB/$PPROB/g" -e "s/CCC/$PPRON/g" TCL/snapshot3.tcl > ___ligb.tcl
sed -e "s/RRR/$FF/g" ___ligb.tcl  > __ligb.tcl
vmd -dispdev text -e __ligb.tcl

mkdir snapshot-$FOUTDIR/zpdb-frame-c$CUTOFF3
mv  pro*.pdb snapshot-$FOUTDIR/zpdb-frame-c$CUTOFF3
mv  lig*.pdb snapshot-$FOUTDIR/zpdb-frame-c$CUTOFF3

done
fi

cat snapshot-$FOUTDIR/zpdb-frame-c$CUTOFF3/lig-$FF-$FFRAM-*.pdb > ____xx
grep -v 'END'    ____xx   > ____xx2
grep -v 'CRYST1' ____xx2  > snapshot-$FOUTDIR/zpdb-frame-c$CUTOFF3/lig-$FF-$FFRAM.pdb

done
rm -r __*
# Extract protein and ligand ##

done
