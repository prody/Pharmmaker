# 2019-06-30 Bahar LAB, University of Pittsburgh  bahar@pitt.edu
# Written by Ji Young Lee
# Calculate interaction statistics from snapshots of 
# Druggability molecular dynamics simulations
# and rank dominant interactions

######### Check here  #######
# PDB file name and path
PDB='struc-list.dat'
# DCD file name and path RRR is variable
DCD='traj-list.dat'
# DCD trajectory number for RRR
TRJNUM=' 1 '
# frame number for each dcd         
FRANUM=10000
# protein chain ID in pdb
CHAIN='chain-list.dat'
# probe molecule ID in pdb
PROBE='probe-list.dat' 
# cutoff between residue and probe      
CUTOFF=4.0
# cutoff between probe and hot spot
CUTOFF2=1.5
# frequency of dcd for analysis
STEP=1
# input directory
INDIR='highaffresid'
# output directory             
OUTDIR=' s1 '
######### Check here  #######

for arg in `seq 1 "$#"`; do
  if [[ "$arg" -eq 1 ]]; then
    PDB=$1
  elif [[ "$arg" -eq 2 ]]; then
    DCD=$2
  elif [[ "$arg" -eq 3 ]]; then
    TRJNUM=$3
  elif [[ "$arg" -eq 4 ]]; then
    FRANUM=$4
  elif [[ "$arg" -eq 5 ]]; then
    CHAIN=$5
  elif [[ "$arg" -eq 6 ]]; then
    PROBE=$6
  elif [[ "$arg" -eq 7 ]]; then
    CUTOFF=$7
  elif [[ "$arg" -eq 8 ]]; then
    CUTOFF2=$8
  elif [[ "$arg" -eq 9 ]]; then
    STEP=$9
  elif [[ "$arg" -eq 10 ]]; then
    INDIR=${10}
  elif [[ "$arg" -eq 11 ]]; then
    OUTDIR=${11}
  fi
done

if [[ $1 = "help" ]]; then
  echo "This script takes up to 7 arguments in the following order"
  echo ""
  echo "protein PDB file path, "
  echo ""
  echo "protein DCD file path, "
  echo ""
  echo "trajectory numbers for multiple runs, "
  echo ""
  echo "number of frames to use, "
  echo ""
  echo "chain ID, "
  echo ""
  echo "probe type, "
  echo ""
  echo "cutoff between residue and probe, "
  echo ""
  echo "cutoff between probe and hot spot, "
  echo ""
  echo "step for reading the trajectories, "
  echo ""
  echo "output directory for results, "
  echo ""
  exit
fi

# Set DCD, PDB, CHAIN and PROBE
if [[ DCD != "*dcd" ]]; then
  DCD=`cat $DCD` || DCD=$DCD
fi
DCD=`echo $DCD | sed 's/ /,/'`

if [[ PDB != "*pdb" ]]; then
  PDB=`cat $PDB` || DCD=$DCD
fi
PDB=`echo $PDB | sed 's/ /,/'`

CHAIN=`cat $CHAIN` || CHAIN=$CHAIN

PROBE=`cat $PROBE` || PROBE=$PROBE

for FOUTDIR in $OUTDIR
do

mkdir -p snapshot-$FOUTDIR

for FCHAIN in $CHAIN
do
for FPROBE in $PROBE
do

grep "$FPROBE $FCHAIN" $INDIR/highaffresid2.dat > ____tt
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
cd __run
env VMDARGS='text with blanks' vmd -dispdev text -e $PHARMMAKER_HOME/snapshot1.tcl -args ../$PDB ../$DCD $STEP $CUTOFF $FPROBE $EE $FCHAIN $FF
cd ..
done

mv  __run/ligbo*  $resdir
cat __run/v*pdb > $resdir/v-com.pdb
ls  __run/v*pdb >  __test
awk '{print $1, (NR-1) }' __test > $resdir/zlist

rm -r __run
###############################

grep $FPROBE $INDIR/hotspots2.pdb > _____test0
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
cd $resdir
#mkdir $kesdir
env VMDARGS='text with blanks' vmd -dispdev text -e $PHARMMAKER_HOME/snapshot2.tcl -args $STEP $CUTOFF2 $FPROBE $FF
#mv ligbo-ok.dat xok*pdb    $kesdir
mv _ligbo-ok.dat ligbo-ok.$FF.dat

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
