#!/bin/bash

CL=8                               # cutoff to include hot spots
hotspotsFile='site-list.dat'       # file/directory path for hotspots PDB
highaffresidsFile='site-list.dat'   # file/directory path for highaffresids.dat
proteinFile='struc-list.dat'       # file/directory for protein PDB
CHAIN='chain-list.dat'             # chain IDs
PROBE='probe-list.dat'             # probe types

if [[ $1 = "help" ]]; then
  echo "This script takes 6 arguments in the following order:"
  echo ""
  echo "path to hotspots file or directory containing it, "
  echo ""
  echo "cutoff length for hot spots being near high affinity residues, "
  echo ""
  echo "path to highaffresid file or directory containing it, "
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
  echo ""
  exit
fi


for arg in `seq 1 "$#"`; do
  if [[ "$arg" -eq 1 ]]; then
    hotspotsFile=$1
  elif [[ "$arg" -eq 2 ]]; then
    CL=$2
  elif [[ "$arg" -eq 3 ]]; then
    highaffresidsFile=$3
  elif [[ "$arg" -eq 4 ]]; then
    proteinFile=$4
  elif [[ "$arg" -eq 5 ]]; then
    CHAIN=$5
  elif [[ "$arg" -eq 6 ]]; then
    PROBE=$6
  fi
done

# Check where site-list.dat exists

[ -f site-list.dat ] && siteListExists=1 || siteListExists=0

# Make sure we have valid hot spots files
#############################################################
if [[ $hotspotsFile = "site-list.dat" ]]; then
  if [[ $siteListExists -eq 1 ]]; then
    hotspotsFile=`cat $hotspotsFile`
  else
    hotspotsFile=''
  fi
elif [[ $hotspotsFile = "all" ]]; then
  hotspotsFile=`du -d 2 | awk '{ print $2 }' | head -n -2`
elif [[ $hotspotsFile = "best" ]]; then
  hotspotsFile='highaffresid/best-site'
fi

if [[ $hotspotsFile = '' ]]; then
  echo "Please provide a path to hot spots file or a site directory to find it in"
  echo ""
  exit
fi

noPDBs=1
IFS=','
numHotspotsEntries=0
hotspotsFiles=()
for file in $hotspotsFile; do
  if [[ -d $file ]]; then
    IFS=$' \t\n'
    hotspotsFilesList=`ls $file*pdb`
    for file2 in $hotspotsFilesList; do
      noPDBs=0
      hotspotsFiles+=($file2)
      ((numHotspotsEntries++))
    done
  else
    if [[ $file == *pdb ]]; then
      noPDBs=0
      hotspotsFiles+=($file)
      ((numHotspotsEntries++))
    fi
  fi
done
hotspotsFileList=${hotspotsFiles[@]} # This gives us back a list

if [[ ${#hotspotsFileList[@]} -eq 0 ]]; then 
  echo "Please provide a path to hot spots file or a valid site directory to find it in"
  echo ""
  exit
fi

if [[ $noPDBs -eq 1 ]]; then
  echo "Hot spots files should end in .pdb"
  echo ""
  exit
fi

# Tell the user what the hot spots file is

IFS=$' \t\n'
echo "hotspotsFile is set to contain"
for file in $hotspotsFileList; do
  echo $file
done
echo ""
#############################################################

# Report CL value

echo "CL is set to $CL"
echo ""

# Sort out highaffresids files
#############################################################
selectedBest=0
if [[ $highaffresidsFile = "site-list.dat" ]]; then
  if [[ $siteListExists -eq 1 ]]; then
    highaffresidsFile=`cat $highaffresidsFile`
  else
    highaffresidsFile='highaffresid'
  fi
elif [[ $highaffresidsFile = "all" ]]; then
  highaffresidsFile=`du -d 2 | awk '{ print $2 }' | head -n -2`
elif [[ $highaffresidsFile = "best" ]]; then
  highaffresidsFile='highaffresid/best-site'
  selectedBest=1
fi

echo "highaffresidsFile is set to contain"
for file in $highaffresidsFile; do
  echo $file
done
echo ""

numHighaffresidDirs=0
highaffresidsFiles=()
for file in $highaffresidsFile; do
  if [[ -d $file ]]; then
    highaffresidsFiles+=(`ls $file/*highaffresid.dat`)
    ((numHighaffresidDirs++))
  else
    highaffresidsFiles+=($file)
  fi
done
highaffresidsFileList=${highaffresidsFiles[@]} # This gives us back a lists

if [[ $numHighaffresidDirs -eq 0 ]]; then
  highaffresidsDirs='highaffresid'
else
  highaffresidsDirs=$highaffresidsFile
fi

# Make sure it gives us highaffresid.dat files

for file in $highaffresidsFileList; do
  if [[ $file != *highaffresid.dat ]]; then
    echo "All highaffresid files should end in highaffresid.dat"
    echo ""
    exit
  fi
done

echo "highaffresidsFileList is set to contain"
for file in $highaffresidsFileList; do
  echo $file
done
echo ""
#############################################################

# Set proteinFile, CHAIN and PROBE

if [[ $proteinFile != "*pdb" ]]; then
  proteinFile=`cat $proteinFile` || proteinFile=$proteinFile
fi

echo "proteinFile is set to contain"
for file in $proteinFile; do
  echo $file
done
echo ""

if [[ -f $CHAIN ]]; then
  CHAIN=`cat $CHAIN` 
else
  IFS=','
  for chain in $CHAIN; do
    chains+=($chain)
  done
  CHAIN=${chains[@]}
fi
IFS=$' \t\n'

echo "CHAIN is set to contain"
for CC in $CHAIN; do
  echo $CC
done
echo ""

if [[ -f $PROBE ]]; then
  PROBE=`cat $PROBE` 
else
  IFS=','
  for probe in $PROBE; do
    probes+=($probe)
  done
  PROBE=${probes[@]}
fi
IFS=$' \t\n'

echo "PROBE is set to contain"
for probe in $PROBE; do
  echo $probe
done
echo ""


echo "There are $numHotspotsEntries hot spots file entries"
echo "There are $numHighaffresidDirs highaffresid file directories"
echo ""

# Start actual calculations

for i in `seq 0 $(expr $numHighaffresidDirs - 1)`; do
  arr=($highaffresidsDirs)
  dir=(${arr[${i}]})
  if [[ $dir = *best-site ]]; then
    if [[ $selectedBest -eq 0 ]]; then
      ((i++))
      continue
    fi
  fi

  for hotspotsFile in $hotspotsFileList; do

    suffix=`echo $hotspotsFile | awk -F/ '{ print $NF }' | awk -F. '{ print $1 }'`
    echo $suffix >> site-list2.dat
    dir2="highaffresid/$suffix"
    mkdir -p $dir2

    for CC in $CHAIN; do
      for probe in $PROBE; do

        for file in $highaffresidsFileList; do
          highaffresidsfile=`echo $file | grep $dir | grep $probe | grep "/$CC"`
          if [ ! -z "$highaffresidsfile" ]; then
            highaffresids=`head $highaffresidsfile`
            if [ ! -z "$highaffresids" ]; then
              echo $highaffresids

              echo "running command: hotspotsnearhighaffresids.sh $CC $probe $CL $highaffresidsfile $proteinFile $hotspotsFile $dir2"
              hotspotsnearhighaffresids.sh $CC $probe $CL $highaffresidsfile $proteinFile $hotspotsFile $dir2
              
              if [[ -f $dir2/$probe-$CC/$probe-$CC-all2-highaffresid.dat ]]; then
                cat $dir2/$probe-$CC/$probe-$CC-all2-highaffresid.dat >> $dir2/highaffresid2.dat
              fi

              if [[ -f $dir2/$probe-$CC/$probe-$CC-all-hs.pdb ]]; then
                cat $dir2/$probe-$CC/$probe-$CC-all-hs.pdb >> $dir2/hotspots.pdb
              fi

              rm -rf $dir2/$probe-$CC/

            else
              echo "No high affinity residues for this condition"
            fi
          fi
        done # end loop for highaffresid files
      done # end loop for probes
    done # end loop for chains

    if [[ -f $dir2/hotspots.pdb ]]; then
      sort --unique $dir2/hotspots.pdb > $dir2/hotspots2.pdb
      sort --unique $dir2/hotspots.pdb > $dir/all-hotspots2.pdb
      rm $dir2/hotspots.pdb
    fi

  done # end loop for hotspots files

done # end loop for directories

echo ""

exit
