# 2019-06-30 Bahar LAB, University of Pittsburgh  bahar@pitt.edu
# Written by Ji Young Lee
# Calculate interaction statistics from snapshots of 
# Druggability molecular dynamics simulations
# and rank dominant interactions

######### Check here  #######
# input directory
INDIR='site-list2.dat'
# cutoff between residue and probe      
CUTOFF=4.0
# cutoff between probe and hot spot
CUTOFF2=1.5
# PDB file name and path
PDB='struc-list.dat'
# DCD file name and path
DCD='traj-list.dat'
# protein chain ID in pdb
CHAIN='chain-list.dat'
# probe molecule ID in pdb
PROBE='probe-list.dat' 
# first frame number for each trajectory 
frameFirst='first'
# last frame number for each trajectory        
frameLast='last'
######### Check here  #######

for arg in `seq 1 "$#"`; do
  if [[ "$arg" -eq 1 ]]; then
    INDIR=$1
  elif [[ "$arg" -eq 2 ]]; then
    CUTOFF=$2
  elif [[ "$arg" -eq 3 ]]; then
    CUTOFF2=$3
  elif [[ "$arg" -eq 4 ]]; then
    PDB=$4
  elif [[ "$arg" -eq 5 ]]; then
    DCD=$5
  elif [[ "$arg" -eq 7 ]]; then
    CHAIN=$7
  elif [[ "$arg" -eq 8 ]]; then
    PROBE=$8
  elif [[ "$arg" -eq 9 ]]; then
    frameFirst=$9
  elif [[ "$arg" -eq 10 ]]; then
    frameLast=${10}
  fi
done

if [[ $1 = "help" ]]; then
  echo "This script takes up to 10 arguments in the following order"
  echo ""
  echo "input directory for results from site selection step"
  echo "The last part of this is used as a suffix for creating a new directory"
  echo ""
  echo "cutoff between residue and probe, "
  echo ""
  echo "cutoff between probe and hot spot, "
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

# Set DCD, PDB, CHAIN and PROBE
if [[ DCD != "*dcd" ]]; then
  DCD=`cat $DCD` || DCD=$DCD
fi
DCD=`echo $DCD | sed 's/ /,/'`

if [[ PDB != "*pdb" ]]; then
  PDB=`cat $PDB` || PDB=$PDB
fi
PDB=`echo $PDB | sed 's/ /,/'`

CHAIN=`cat $CHAIN` || CHAIN=$CHAIN

PROBE=`cat $PROBE` || PROBE=$PROBE

if [[ -f $INDIR ]]; then
  INDIR=`cat $INDIR`
fi

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
  mkdir -p $FOUTDIR

  for FCHAIN in $CHAIN
  do
    for FPROBE in $PROBE
    do

    grep "$FPROBE $FCHAIN" highaffresid/$FINDIR/highaffresid2.dat > ____tt
    sed -e "s/$FPROBE $FCHAIN/    /g" ____tt > ____tt1
    RESID=`cat ____tt1 | awk '{print $0}'`


    for EE in $RESID
    do
      resdir=z.$FCHAIN.$EE.$FPROBE
      mkdir $resdir

      mkdir __run
      cd __run
      env VMDARGS='text with blanks' vmd -dispdev text -e $PHARMMAKER_HOME/snapshot1.tcl -args ../$PDB ../$DCD $CUTOFF $FPROBE $EE $FCHAIN
      cd ..

      mv  __run/ligbo*  $resdir/
      cat __run/v*pdb > $resdir/v-com.pdb
      ls  __run/v*pdb >  __test
      awk '{print $1, (NR-1) }' __test > $resdir/zlist

      rm -r __run
      ###############################

      grep $FPROBE highaffresid/$FINDIR/hotspots2.pdb | sort -n -k10 > __test1

      TNUM=`wc -l __test1 | awk '{print $1}'`

      echo "/END/{i\\" > __test3

      for (( mm = 1 ; mm <= $TNUM  ; mm++ ))
      do
        sed -n -e "$mm,$mm p" __test1 > __test2

        awk '{printf("%-6s%5s  %-3s %4s%1s%4d %11.3f%8.3f%8.3f%6.2f%6.2f%13s\n", "ATOM", "'"'$mm'"'", "C2", "'"'$FPROBE'"'", "M", "'"'$mm'"'", $6, $7, $8, $9, $10, "M\\" )}' \
        __test2 >> __test3
        awk '{printf("%-6s%5s  %-3s %4s%1s%4d %11.3f%8.3f%8.3f%6.2f%6.2f%13s\n", "ATOM", "'"'$mm'"'", "C2", "'"'$FPROBE'"'", "M", "'"'$mm'"'", $6, $7, $8, $9, $10, "M" )}' \
        __test2 >> __ttest3
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

        env VMDARGS='text with blanks' vmd -dispdev text -e $PHARMMAKER_HOME/snapshot2.tcl -args $CUTOFF2 $FPROBE $FF

        mv _ligbo-ok.dat ligbo-ok.$FF.dat

        TNNN=`wc -l zlist | awk '{print $1}'`

        for (( rr = 0   ; rr <= $TNNN     ; rr++ ))
        do

          awk '{if ($1 == "'"$rr"'") print }' ligbo-ok.$FF.dat > ___test2

          CRAT=`wc -l ___test2 | awk '{print $1}'`

          if [ $CRAT -ge 1 ]; then
            awk '{if ($2 == "'"$rr"'") print }' zlist  > ___test3

            perl -i -pe 's/-/ /g' ___test3
            perl -i -pe 's/.pdb/    /g' ___test3

            FRAMT=`cat ___test3 | awk '{print $2}'`
            FRAME=`cat ___test3 | awk '{print $3}'`
            FRAMM=`cat ___test3 | awk '{print $4}'`

            awk '{if ($1 == "'"$FRAMM"'") print "'"$FRAMT"'" , "'"$FRAME"'" , $3, $4, $5, $7, $8, $10 }' ligbo-ok.$FF.dat >> tout-hotspot-$FF.dat
            awk '{if ($1 == "'"$FRAME"'") print "'"$FRAMT"'" , $0 }' ligbo-$FRAMT.dat        >> tout-res-$FF.dat

            awk '{print "SIM#", $2, "FRAME#", $3 }' ___test3 >> tout-$FF.dat
          fi

        done

        rm _*
        cd ..

        #### Sort
        KK=0
        for traj in $DCD; do
          awk '{if ($2 == "'"$KK"'") print }' $resdir/tout-$FF.dat | sort -n -k4 >> $resdir/out-$FF.dat
          awk '{if ($1 == "'"$KK"'") print }' $resdir/tout-res-$FF.dat | sort -n -k2 >> $resdir/out-res-$FF.dat
          awk '{if ($1 == "'"$KK"'") print }' $resdir/tout-hotspot-$FF.dat | sort -n -k2 >> $resdir/out-hotspot-$FF.dat 
          ((KK++))
        done
        #### Sort

        wc -l $resdir/out-$FF.dat  >> $resdir/out-count
        cat $resdir/out-$FF.dat  >> ___tttt            

      done

      #########

      rm    $resdir/v*pdb $resdir/lig* $resdir/z* $resdir/t*

      ######### outfr.dat
      if [[ "$frameFirst" -eq "first" ]]; then
        frameFirst=1
      fi

      if [[ "$frameLast" -eq "last" ]]; then
        frameLast=`tail -n1 $resdir/out-?.dat | sort -n -k4 | tail -n1 | awk '{ print $NF }'`
      fi

      KK=0
      for traj in $DCD
      do
        for (( nn = $frameFirst ; nn <= $frameLast   ; nn++ ))
        do
          CRAT=`awk '{if ($2 == "'"$KK"'" && $4 == "'"$nn"'" ) print }' ___tttt | wc -l | awk '{print $1}'`

          if [ $CRAT -ge 1 ]; then
            PROBID=`grep "$KK $nn" $resdir/out-res-*.dat | head -n 1 | awk '{print $8}'`
            PROBNU=`grep "$KK $nn" $resdir/out-res-*.dat | head -n 1 | awk '{print $9}'`
            echo "SIM# $KK FRAME# $nn PROBE $PROBID PROBE# $PROBNU "  >> $resdir/outfr.dat
          fi
        done
        ((KK++))
      done
      rm _*
      #########

      mv $resdir $FOUTDIR/

      done
      rm _*
    done

  done

  wc -l $FOUTDIR/z.*.*/outfr* | grep -v "total" | sort -r -n -k1 > $FOUTDIR/zlist-count

  rm __*

done
