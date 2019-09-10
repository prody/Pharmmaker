# 2019-06-30 Bahar LAB, University of Pittsburgh  bahar@pitt.edu
# Written by Ji Young Lee
# Calculate interaction statistics from snapshots of 
# Druggability molecular dynamics simulations
# and rank dominant interactions

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
CHAIN=' A B '
# probe molecule ID in pdb
PROBE=' IPRO ACAM ACTT IPAM IBTN IMID '
# cutoff between residue and probe      
CUTOFF=4.0
# cutoff between probe and hot spot
CUTOFF2=1.5
# frequency of dcd for analysis
STEP=1
# output directory             
OUTDIR=' s1 '
######### Check here  #######

for FOUTDIR in $OUTDIR
do
for FCHAIN in $CHAIN
do
for FPROBE in $PROBE
do

grep "$FPROBE $FCHAIN" ./snapshot-$FOUTDIR/input-highaffresid > ____tt
sed -e "s/$FPROBE $FCHAIN/    /g" ____tt > ____tt1
RESID=`cat ____tt1 | awk '{print $0}'`

# Check Current EE 154 = Tyr151
for EE in $RESID
do
resdir=z.$FCHAIN.$EE.$FPROBE
mkdir $resdir

# Check Current FF ############
mkdir __run
for FF in $TRJNUM
do
sed -e "s/APDB/$PDB/g" -e "s/ADCD/$DCD/g" -e "s/SSTEP/$STEP/g" -e "s/CUTOFF/$CUTOFF/g" -e "s/AAA/$FPROBE/g" -e "s/BBB/$EE/g" -e "s/CCC/$FCHAIN/g" TCL/snapshot1.tcl > ___ligb.tcl      
sed -e "s/RRR/$FF/g" ___ligb.tcl  > __ligb.tcl
cd __run
vmd -dispdev text -e ../__ligb.tcl
cd ..
done

mv  __run/ligbo*  $resdir
cat __run/v*pdb > $resdir/v-com.pdb
ls  __run/v*pdb >  __test
awk '{print $1, (NR-1) }' __test > $resdir/zlist

rm -r __run
###############################

grep $FPROBE ./snapshot-$FOUTDIR/hotspot.pdb > _____test0
sort -n -k10 _____test0             > __test1

TNUM=`wc -l __test1 | awk '{print $1}'`

echo "/END/{i\\" > __test3

for (( mm = 1 ; mm <= $TNUM  ; mm++ ))
do
sed -n -e "$mm,$mm p" __test1 > __test2

cat > ___molexe  <<EOF
awk '{printf("%-6s%5s  %-3s %4s%1s%4d %11.3f%8.3f%8.3f%6.2f%6.2f%13s\n", "ATOM", "$mm", "C2", "$FPROBE", "M", "$mm", \$6, \$7, \$8, \$9, \$10, "M\\\" )}' \
__test2 >> __test3
awk '{printf("%-6s%5s  %-3s %4s%1s%4d %11.3f%8.3f%8.3f%6.2f%6.2f%13s\n", "ATOM", "$mm", "C2", "$FPROBE", "M", "$mm", \$6, \$7, \$8, \$9, \$10, "M" )}' \
__test2 >> __ttest3
EOF
chmod 777 ___molexe
./___molexe
rm ___molexe
done

cp __ttest3 $resdir/hot-spot.pdb

echo "END " >> __test3
echo "d "   >> __test3
echo "} "   >> __test3

sed -f __test3  $resdir/v-com.pdb > $resdir/v-com-ok.pdb
rm _*

#########
for (( FF = 1 ; FF <= $TNUM  ; FF++ ))
do
sed -e "s/SSTEP/$STEP/g" -e "s/CUTOFF/$CUTOFF2/g" -e "s/AAA/resname $FPROBE and chain P and not hydrogen/g"  -e "s/BBB/resname $FPROBE and chain M and resid $FF/g" TCL/snapshot2.tcl > $resdir/__ligb.tcl
cd $resdir
#mkdir $kesdir
vmd -dispdev text -e __ligb.tcl
#mv ligbo-ok.dat xok*pdb    $kesdir
mv _ligbo-ok.dat ligbo-ok.$FF.dat
rm __ligb.tcl


TNNN=`wc -l zlist | awk '{print $1}'`

for (( rr = 0   ; rr <= $TNNN     ; rr++ ))
do

cat > _molexe  <<EOF
awk '{if (\$1 == $rr) print }' ligbo-ok.$FF.dat > ___test2
EOF
chmod 777 _molexe
./_molexe
rm _molexe

CRAT=`wc -l ___test2 | awk '{print $1}'`

if [ $CRAT -ge 1 ]; then
cat > _molexe  <<EOF
awk '{if (\$2 == $rr) print }' zlist  > ___test3
EOF
chmod 777 _molexe
./_molexe
rm _molexe

perl -i -pe 's/-/ /g' ___test3
perl -i -pe 's/.pdb/    /g' ___test3

FRAMT=`cat ___test3 | awk '{print $2}'`
FRAME=`cat ___test3 | awk '{print $3}'`
FRAMM=`cat ___test3 | awk '{print $4}'`

cat > _molexe  <<EOF
awk '{if (\$1 == $FRAMM) print $FRAMT , $FRAME , \$3, \$4, \$5, \$7, \$8, \$10 }' ligbo-ok.$FF.dat >> tout-hotspot-$FF.dat
awk '{if (\$1 == $FRAME) print $FRAMT , \$0 }' ligbo-$FRAMT.dat        >> tout-res-$FF.dat
EOF
chmod 777 _molexe
./_molexe
rm _molexe

awk '{print "SIM#", $2, "FRAME#", $3 }' ___test3 >> tout-$FF.dat
#cat xok-$rr.pdb  >> x-com.$FF.pdb
fi

done

rm _*
cd ..

#### Sort
for KK in $TRJNUM
do
cat > _molexe  <<EOF
awk '{if (\$2 == $KK) print }' $resdir/tout-$FF.dat         > ______test1
awk '{if (\$1 == $KK) print }' $resdir/tout-res-$FF.dat     > ______test2
awk '{if (\$1 == $KK) print }' $resdir/tout-hotspot-$FF.dat > ______test3
EOF
chmod 777 _molexe
./_molexe
rm _molexe

sort -n -k4 ______test1 >>  $resdir/out-$FF.dat         
sort -n -k2 ______test2 >>  $resdir/out-res-$FF.dat     
sort -n -k2 ______test3 >>  $resdir/out-hotspot-$FF.dat 
done
#### Sort

wc -l $resdir/out-$FF.dat  >> $resdir/out-count
cat $resdir/out-$FF.dat  >> ___tttt            

done

#########
rm    $resdir/v*pdb $resdir/lig* $resdir/z* $resdir/t*

######### outfr.dat
for KK in $TRJNUM
do
for (( nn = 0   ; nn <= $FRANUM   ; nn++ ))
do
cat > _molexe  <<EOF
awk '{if (\$2 == $KK && \$4 == $nn ) print }' ___tttt  > ___tttt2
EOF
chmod 777 _molexe
./_molexe
rm _molexe

CRAT=`wc -l ___tttt2 | awk '{print $1}'`

if [ $CRAT -ge 1 ]; then
grep "$KK $nn" $resdir/out-res-*.dat > ____tttt
PROBID=`head -n 1 ____tttt | awk '{print $8}'`
PROBNU=`head -n 1 ____tttt | awk '{print $9}'`
echo "SIM# $KK FRAME# $nn PROBE $PROBID PROBE# $PROBNU "  >> ___tttt3
fi
done
done
mv ___tttt3 $resdir/outfr.dat
rm _*
#########

mv $resdir snapshot-$FOUTDIR

done
rm _*
done
done

wc -l snapshot-$FOUTDIR/z.*.*/outfr*  >  ___tttt4
grep -v "total" ___tttt4              >  ___tttt5
sort -r -n -k1  ___tttt5              >  ___tttt6
cp ___tttt6 snapshot-$FOUTDIR/zlist-count

rm __*

done
