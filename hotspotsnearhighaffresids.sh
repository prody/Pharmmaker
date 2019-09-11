# 2015-05-22 Ji Young Lee
#!/bin/bash

CL=8                # R: redius to include hot spots
CC=B                # chain ID
probe=IBTN         # output directory
highaffresidsfile="../highaffresid/out-$CC-$probe-highaffresid.dat"
proteinFile="../dg/dg_protein_heavyatoms.pdb"
hotspotsFile="../dg/dg_all_hotspots.pdb"
dir="."

for arg in `seq 1 "$#"`; do
  if [[ "$arg" -eq 1 ]]; then
    CC=$1
  elif [[ "$arg" -eq 2 ]]; then
    probe=$2
  elif [[ "$arg" -eq 3 ]]; then
    CL=$3
  elif [[ "$arg" -eq 4 ]]; then
    highaffresidsfile=$4
  elif [[ "$arg" -eq 5 ]]; then
    proteinFile=$5
  elif [[ "$arg" -eq 6 ]]; then
    hotspotsFile=$6
  elif [[ "$arg" -eq 7 ]]; then
    dir=$7
  fi
done

if [[ $1 = "help" ]]; then
  echo "This script takes up to 7 arguments in the following order"
  echo "chain ID, probe type, cutoff for including probes in hotspots, "
  echo "highaffresid file path, protein PDB file path, "
  echo "hotspots PDB file path, and directory for putting the results"
  exit
fi

highaffresids=`cat $highaffresidsfile`

[ ! -d "$dir/$probe-$CC" ] && mkdir $dir/$probe-$CC

### check  residues for probe

grep -v 'REMARK' $hotspotsFile | grep "$probe" > __test2
TNUM2=`wc -l __test2 | awk '{print $1}'`

for MM in $highaffresids; do
  NN=$(printf %4s $MM)
  grep " $CC " $proteinFile | awk -F "" '{ if ( $23 $24 $25 $26 == "'"$NN"'" )  print }' > __test1
  TNUM1=`wc -l __test1 | awk '{print $1}'`

  ### check hot spots trajectory
  for KK in  1 # "1-6"
  do

    # Loop through lines for that resid
    for (( mm = 1  ; mm <= $TNUM1 ; mm++ )); do

      # get coordinates
      X1=`sed -n "${mm},${mm}p" __test1 | awk -F "" '{print $31 $32 $33 $34 $35 $36 $37 $38 }'`
      Y1=`sed -n "${mm},${mm}p" __test1 | awk -F "" '{print $39 $40 $41 $42 $43 $44 $45 $46 }'`
      Z1=`sed -n "${mm},${mm}p" __test1 | awk -F "" '{print $47 $48 $49 $50 $51 $52 $53 $54 }'`
      
      # Loop through lines for hotspots of that probe type
      for (( pp = 1 ; pp <= $TNUM2 ; pp++ )); do

        # get coordinates
        X2=`sed -n "${pp},${pp}p" __test2 | awk -F "" '{print $31 $32 $33 $34 $35 $36 $37 $38 }'`
        Y2=`sed -n "${pp},${pp}p" __test2 | awk -F "" '{print $39 $40 $41 $42 $43 $44 $45 $46 }'`
        Z2=`sed -n "${pp},${pp}p" __test2 | awk -F "" '{print $47 $48 $49 $50 $51 $52 $53 $54 }'`
        
        # calculate distance as float with bc
        xd2=`echo "($X1-$X2)*($X1-$X2)" | bc`
        yd2=`echo "($Y1-$Y2)*($Y1-$Y2)" | bc`
        zd2=`echo "($Z1-$Z2)*($Z1-$Z2)" | bc`
        d2=`echo "($xd2 + $yd2 + $zd2)" | bc`
        distance=`echo "sqrt($d2)" | bc`

        # use bc for float comparison
        st=`echo "$distance <= $CL" | bc`
        if [[ $st -eq 1 ]]; then

          # write pairs to dat
          echo "`sed -n "${mm},${mm}p" __test1 | awk -F "" '{ for(i=12; i<=26; ++i) printf "%s", $i }'` `sed -n "${pp},${pp}p" __test2 | awk -F "" '{ for(i=1; i<=26; ++i) printf "%s", $i }'`     $distance" >> $dir/hotspot-highaffresid-pairs.dat
          
          # initialise highaffresid.dat with probe and chain string
          if [[ ! -f $dir/$probe-$CC/$probe-$CC-all2-highaffresid.dat ]]; then
            printf "$probe $CC " > $dir/$probe-$CC/$probe-$CC-all2-highaffresid.dat
          fi

          # write intermediate files for 
          sed -n "${pp},${pp}p" __test2 >> $dir/$probe-$CC/out-$MM-t$KK-hs.pdb
          sed -n "${mm},${mm}p" __test1 | awk '{ print $6 }' >> $dir/$probe-$CC/out-$MM-t$KK-highaffresid.dat

        fi

      done # finish loop of hotspots

      # extract unique hotspots to the main file
      if [[ -f $dir/$probe-$CC/out-$MM-t$KK-hs.pdb ]]; then
        awk '!seen[$0]++' $dir/$probe-$CC/out-$MM-t$KK-hs.pdb >> $dir/$probe-$CC/$probe-$CC-all-hs.pdb
      fi

    done # finish loop over highaffresid atoms

  done # finish loop of files

  # extract unique highaffresids to main file with formatting
  if [[ -f $dir/$probe-$CC/out-$MM-t$KK-highaffresid.dat ]]; then
    dt=`awk '!seen[$0]++' $dir/$probe-$CC/out-$MM-t$KK-highaffresid.dat` 
    dt=${dt/ /$'\n'/}
    echo -e "$dt \c"  >> $dir/$probe-$CC/$probe-$CC-all2-highaffresid.dat
  fi

done # finish loop over highaffresids

# Make a newline at the end of the line for that that probe type
if [[ -f $dir/$probe-$CC/$probe-$CC-all2-highaffresid.dat ]]; then
  echo "" >> $dir/$probe-$CC/$probe-$CC-all2-highaffresid.dat
fi

rm __*
