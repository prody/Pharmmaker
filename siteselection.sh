#!/bin/bash

CL=8                               # cutoff to include hot spots
hotspotsFile='all'                 # file/directory path for hotspots PDB
highaffresidsFile='highaffresid'   # file/directory path for highaffresids.dat
proteinFile='struc-list.dat'       # file/directory for protein PDB
CHAIN='chain-list.dat'             # chain IDs
PROBE='probe-list.dat'             # probe types
#dir is not set as it will follow highaffresidsFile if the user doesn't provide it

if [[ $1 = "help" ]]; then
  echo "This script takes 6 arguments in the following order:"
  echo ""
  echo "cutoff length for hot spots being near high affinity residues, "
  echo ""
  echo "path to hotspots file or directory containing it, "
  echo ""
  echo "path to highaffresid file or directory containing it, "
  echo ""
  echo "path to directory for output, "
  echo ""
  echo "path to PDB structure or directory containing it, "
  echo ""
  echo "list of chains to be used, "
  echo ""
  echo "and list of probes to be used"
  echo ""
  echo "These values can all be provided inside files such as chain-list.dat"
  echo ""
  echo "The default behaviour is to use the files from the previous step."
  echo "When results are calculated for multiple sites/regions, we use the best one by default."
  echo "All of them can be used by providing the word 'all' for hotspots file."
  echo "Please note that if no value is provided for highaffresid path, the hotspots path value is used."
  exit
fi


for arg in `seq 1 "$#"`; do
  if [[ "$arg" -eq 1 ]]; then
    CL=$1
  elif [[ "$arg" -eq 2 ]]; then
    hotspotsFile=$2
  elif [[ "$arg" -eq 3 ]]; then
    highaffresidsFile=$3
  elif [[ "$arg" -eq 4 ]]; then
    dir=$4
  elif [[ "$arg" -eq 5 ]]; then
    proteinFile=$5
  elif [[ "$arg" -eq 6 ]]; then
    CHAIN=$6
  elif [[ "$arg" -eq 7 ]]; then
    PROBE=$7
  fi
done

echo "CL is set to $CL"
echo ""

# Set hotspotsFile to point to actual directory

echo "hotspotsFile is set to $hotspotsFile"

if [[ $hotspotsFile = "all" ]]; then
  hotspotsFile=`du -d 1 | awk '{ print $2 }' | head -n -1`
fi

if [[ $hotspotsFile = "best" ]]; then
  hotspotsFile='best-site'
fi

echo "hotspotsFile is set to contain"
for file in $hotspotsFile; do
  echo $file
done
echo ""

# Set highaffresids to point to actual files

if [ -z ${highaffresidsFile+x} ]; then 
  highaffresidsFile=$hotspotsFile
else
  if [[ $highaffresidsFile = "all" ]]; then
    highaffresidsFile=`du -d 1 | awk '{ print $2 }' | head -n -1`
  fi

  if [[ $highaffresidsFile = "best" ]]; then
    highaffresidsFile='best-site'
  fi
fi

highaffresidsDirs=()
isDir=0
for file in $highaffresidsFile; do
  if [[ -d $file ]]; then
    highaffresidsDirs+=($file)
    isDir=1
  else
    break
  fi
done

if [[ "$isDir" -eq 1 ]]; then
  highaffresidsFileList=${highaffresidsDirs[@]} # This gives us back a list
else
  highaffresidsFileList=$highaffresidsFile
fi

echo "highaffresidsFileList is set to contain"
for file in $highaffresidsFileList; do
  echo $file
done

if [[ ${#highaffresidsDirs[@]} -eq 0 ]]; then
  highaffresidsDirs+=$(dirname "$highaffresidsFile")
  num_haf_dirs=1
else
  num_haf_dirs=${#highaffresidsDirs[@]}
fi

echo "the array contains $num_haf_dirs directories"
echo ""

highaffresidsFiles=()
for file in $highaffresidsFile; do
  if [[ -d $file ]]; then
    highaffresidsFiles+=(`ls $file/*highaffresid.dat`)
  else
    highaffresidsFiles+=($file)
  fi
done
highaffresidsFileList=${highaffresidsFiles[@]} # This gives us back a list

echo "highaffresidsFileList is set to contain"
for file in $highaffresidsFileList; do
  echo $file
done
echo ""

# Set hotspotsFile to point to actual files (after passing on the directories)

hotspotsFiles=()
for file in $hotspotsFile; do
  if [[ -d $file ]]; then
    hotspotsFiles+=(`ls $file/*site_?_soln_?.pdb`)
  else
    hotspotsFiles+=($file)
  fi
done
hotspotsFileList=${hotspotsFiles[@]} # This gives us back a list

echo "hotspotsFileList is set to contain"
for file in $hotspotsFileList; do
  echo $file
done
num_hs_files=${#hotspotsFiles[@]}
echo "the array contains $num_hs_files filenames"
echo ""

#if [[ $num_hs_files != $num_haf_dirs ]]; then
#  echo "ERROR: the numbers of hotspot files ($num_hs_files) and highaffresid directories ($num_haf_dirs) do not match"
#  echo ""
#  exit
#fi

# Set proteinFile, CHAIN and PROBE

if [[ proteinFile != "*pdb" ]]; then
  proteinFile=`cat $proteinFile` || proteinFile=$proteinFile
fi

echo "proteinFile is set to contain"
for file in $proteinFile; do
  echo $file
done
echo ""

CHAIN=`cat $CHAIN` || CHAIN=$CHAIN
echo "CHAIN is set to contain"
for CC in $CHAIN; do
  echo $CC
done
echo ""

PROBE=`cat $PROBE` || PROBE=$PROBE
echo "PROBE is set to contain"
for probe in $PROBE; do
  echo $probe
done
echo ""

# report about setting dir

if [ -z ${dir+x} ]; then 
  echo "dir is unset and will be set in the loop"
else 
  echo "dir is set to '$dir'"
fi


for i in `seq 0 $(expr $num_haf_dirs - 1)`; do


  if [ -z ${dir+x} ]; then
    dir=${highaffresidsDirs[${i}]}
  fi
  
  j=0
  for hotspotsFile in $hotspotsFileList; do
    ((j++))
    if [ "$j" -gt "$i" ]; then
      break
    fi
  done

  for CC in $CHAIN; do
    for probe in $PROBE; do

      echo ""
      echo "iteration $i, directory $dir, chain $CC, probe $probe"

      for file in $highaffresidsFileList; do
        highaffresidsfile=`echo $file | grep $dir | grep $probe | grep "/$CC"`
        if [ ! -z "$highaffresidsfile" ]; then
          highaffresids=`head $highaffresidsfile`
          if [ ! -z "$highaffresids" ]; then
            echo $highaffresids

            echo "running command: hotspotsnearhighaffresids.sh $CC $probe $CL $highaffresidsfile $proteinFile $hotspotsFile $dir"
            hotspotsnearhighaffresids.sh $CC $probe $CL $highaffresidsfile $proteinFile $hotspotsFile $dir
            
            if [[ -f $dir/$probe-$CC/$probe-$CC-all2-highaffresid.dat ]]; then
              cat $dir/$probe-$CC/$probe-$CC-all2-highaffresid.dat >> $dir/highaffresid2.dat
            fi

            if [[ -f $dir/$probe-$CC/$probe-$CC-all-hs.pdb ]]; then
              cat $dir/$probe-$CC/$probe-$CC-all-hs.pdb >> hotspots.pdb
            fi

            rm -rf $dir/$probe-$CC/

          else
            echo "No high affinity residues for this condition"
          fi
        fi
      done # end loop for highaffresid files
    done # end loop for probes
  done # end loop for chains
done # end loop for directories

sort --unique hotspots.pdb > $dir/hotspots2.pdb
rm hotspots.pdb

echo ""

exit
