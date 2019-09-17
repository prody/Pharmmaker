# 2019-06-30 Bahar LAB, University of Pittsburgh  bahar@pitt.edu
# Written by Ji Young Lee
# Select snapshots with many of the most dominant interactions from
# Druggability molecular dynamics simulations

######### Check here  #######
# Frequency cutoff for selecting interactions
CUTOFF3=0.1
# input directory
INDIR='site-list2.dat'
# PDB file name and path
PDB='struc-list.dat'
# DCD file name and path
DCD='traj-list.dat'
# first frame number for each dcd 
frameFirst='first-frame.dat'
# last frame number for each dcd        
frameLast='last-frame.dat'
######### Check here  #######

if [[ $# -lt 1 ]]; then
  echo "Please provide a frequency score cutoff"
  echo ""
  exit
fi

for arg in `seq 1 "$#"`; do
  if [[ "$arg" -eq 1 ]]; then
    CUTOFF3=$1
  elif [[ "$arg" -eq 2 ]]; then
    INDIR=$2
  elif [[ "$arg" -eq 3 ]]; then
    PDB=$3
  elif [[ "$arg" -eq 4 ]]; then
    DCD=$4
  elif [[ "$arg" -eq 5 ]]; then
    CHAIN=$5
  elif [[ "$arg" -eq 6 ]]; then
    PROBE=$6
  elif [[ "$arg" -eq 7 ]]; then
    frameFirst=$7
  elif [[ "$arg" -eq 8 ]]; then
    frameLast=$8
  fi
done

if [[ $1 = "help" ]]; then
  echo "This script takes up to 8 arguments in the following order"
  echo ""
  echo "cutoff frequency for interactions for selecting snapshots, "
  echo ""
  echo "input directory for results from site selection step "
  echo "This is used to locate the results from the snapshot statistics step"
  echo ""
  echo "protein PDB file path, "
  echo ""
  echo "protein DCD file path, "
  echo ""
  echo "chain ID, "
  echo ""
  echo "probe type, "
  echo ""
  echo "first frame to use, "
  echo ""
  echo "last frame to use, "
  echo ""
  exit
fi

# Set DCD, PDB, INDIR, CHAIN and PROBE
if [[ DCD != "*dcd" ]]; then
  DCD=`cat $DCD` || DCD=$DCD
fi

if [[ PDB != "*pdb" ]]; then
  PDB=`cat $PDB` || DCD=$DCD
fi

if [[ -f $INDIR ]]; then
  INDIR=`cat $INDIR`
fi

if [[ -f $frameFirst ]]; then
  frameFirst=`cat $frameFirst` 
else
  IFS=','
fi

firstFramesArray=()
for firstFrame in $frameFirst; do
  firstFramesArray+=($firstFrame)
done
IFS=$' \t\n'

if [[ -f $frameLast ]]; then
  frameLast=`cat $frameLast` 
else
  IFS=','
fi

lastFramesArray=()
for lastFrame in $frameLast; do
  lastFramesArray+=($lastFrame)
done
IFS=$' \t\n'

IFS=','
INDIRS=()
for indir in $INDIR; do
  INDIRS+=($indir)
done
INDIR=${INDIRS[@]}
IFS=$' \t\n'

for FINDIR in $INDIR
do

  FOUTDIR="snapshot/`echo $FINDIR | awk -F/ '{print $NF}'`"

  #########
  awk -v cutoff=$CUTOFF3 '{if ($1 >= cutoff ) print }' $FOUTDIR/zlist-frequency > ___test 

  mkdir CUT$CUTOFF3

  TNUM=`wc -l ___test | awk '{print $1}'`
  for (( mm = 1 ; mm <= $TNUM  ; mm++ ))
  do
    FF=`sed -n -e "$mm,$mm p" ___test | sed -e "s/outfr.dat/  /g" | awk '{print $2}'`
    cp -r $FF CUT$CUTOFF3
  done

  rm ___t*
  ############

  ######### zlist-frame-c1000 zlist-frame-c1000-detail


  trajNum=1
  for traj in $DCD
  do
    for (( nn = ${firstFramesArray[$trajNum]} ; nn <= ${lastFramesArray[$trajNum]}   ; nn++ ))
    do

      grep "SIM# $trajNum FRAME# $nn "  CUT$CUTOFF3/z.*/outfr.dat > ___test

      CRAT=`wc -l ___test | awk '{print $1}'`

      if [ $CRAT -eq $TNUM ]; then
        cat   ___test >> $FOUTDIR/zlist-frame-c$CUTOFF3-detail
        echo "SIM# $trajNum FRAME# $nn   $CRAT"  >> $FOUTDIR/zlist-frame-c$CUTOFF3
      fi

    done
    ((trajNum++))
  done

  rm  __*
  \rm -r CUT$CUTOFF3
  ######### zlist-frame-c1000 zlist-frame-c1000-detail


  # Extract protein and ligand ##

  SNAPSNUM=`wc -l $FOUTDIR/zlist-frame-c$CUTOFF3 | awk '{print $1}'`

  for (( yy = 1 ; yy <= $SNAPSNUM  ; yy++ ))
  do

    sed -n -e "$yy,$yy p" $FOUTDIR/zlist-frame-c$CUTOFF3  > _____tt
    trajNum=`cat _____tt | awk '{print $2}'`
    FFRAM=`cat _____tt | awk '{print $4}'`
    grep "$trajNum FRAME# $FFRAM" $FOUTDIR/zlist-frame-c$CUTOFF3-detail > _____tt2

    TTTTNUM=`wc -l _____tt2 | awk '{print $1}'`

    if [ $TTTTNUM -ge 1 ]; then
      for (( bb = 1 ; bb <= $TTTTNUM  ; bb++ ))
      do
        sed -n -e  "$bb,$bb p" _____tt2 > _____tt3

        PPROB=`cat _____tt3 | awk '{print $6}'`
        PPRON=`cat _____tt3 | awk '{print $8}'`
        env VMDARGS='text with blanks' vmd -dispdev text -e $PHARMMAKER_HOME/snapshot3.tcl -args $PDB $DCD $FFRAM $PPROB $PPRON $trajNum

        mkdir -p $FOUTDIR/zpdb-frame-c$CUTOFF3
        mv  pro*.pdb $FOUTDIR/zpdb-frame-c$CUTOFF3
        mv  lig*.pdb $FOUTDIR/zpdb-frame-c$CUTOFF3

      done
    fi

    cat $FOUTDIR/zpdb-frame-c$CUTOFF3/lig-$trajNum-$FFRAM-*.pdb > ____xx
    grep -v 'END'    ____xx   > ____xx2
    grep -v 'CRYST1' ____xx2  > $FOUTDIR/zpdb-frame-c$CUTOFF3/lig-$trajNum-$FFRAM.pdb
    rm $FOUTDIR/zpdb-frame-c$CUTOFF3/lig-$trajNum-$FFRAM-*.pdb 

  done
  rm -r __*
  # Extract protein and ligand ##

done
exit
