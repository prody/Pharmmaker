argString=""
i=0
for (( i = 1  ; i <= $# ; i++ )); do
  argString+="${!i} "
done

mkdir -p highaffresid
cd highaffresid
env VMDARGS='text with blanks' vmd -dispdev text -e $PHARMMAKER_HOME/highaffresid.tcl -args $argString
cd ..
