# 2019-06-30 Bahar LAB, University of Pittsburgh  bahar@pitt.edu
# Written by Ji Young Lee
# Select hot spots near high affinity residues from the results of
# DruGUI analyis of druggability molecular dynamics simulations

CL=8                # R: redius to include hot spots
CC=B                # chain ID
outdir=IBTN         # output directory

mkdir $outdir-$CC

### check  residues for probe
for MM in 46 47 105 108 229 232 233 236
do
sed -e "s/-/ -/g"       ../hotspot/dg_protein_heavyatoms.pdb  > __test00
grep " $CC " __test00 > __test0
cat > ___molexe  <<EOF
awk '{if ( \$6 == $MM ) print }' __test0 > __test1
EOF
chmod 777 ___molexe
./___molexe
rm ___*

### check hot spots trajectory
for KK in  1 # "1-6"
do
sed '1,7d'         ../hotspot/dg/dg_all_hotspots.pdb > __test2
sed -e "s/-/ -/g"       __test2     > __test3
cp __test3 __test7
TNUM=`wc -l __test1 | awk '{print $1}'`

####LOOP#######
for (( mm = 1  ; mm <= $TNUM ; mm++ ))
do
sed -n -e "$mm,$mm p" __test1  > __test4
XX=`awk '{print $7}' __test4`
YY=`awk '{print $8}' __test4`
ZZ=`awk '{print $9}' __test4`

cat > ___molexe  <<EOF
awk '{if ( sqrt( (\$6-($XX))*(\$6-($XX)) + (\$7-($YY))*(\$7-($YY)) + (\$8-($ZZ))*(\$8-($ZZ)) ) <= $CL ) print }' __test7 >> __test5
awk '{if ( sqrt( (\$6-($XX))*(\$6-($XX)) + (\$7-($YY))*(\$7-($YY)) + (\$8-($ZZ))*(\$8-($ZZ)) ) > $CL ) print }'  __test7  > __test6
EOF
chmod 777 ___molexe
./___molexe
rm ___*

cp __test6 __test7

done
####LOOP#######

cat > ___molexe  <<EOF
awk '{printf("%-6s%5s  %-3s%6s %3d %11.3f%8.3f%8.3f%6.2f%6.2f           M\n", \$1, \$2, \$3, \$4, \$5, \$6, \$7, \$8, \$9, \$10  )}' \
 __test5 > __test8            
EOF
chmod 777 ___molexe
./___molexe
rm ___*
rm __test5

cp  __test8 ./$outdir-$CC/out-$MM-t$KK.pdb
cat __test8 >> __out-all
cat __test8 >> __out-all2

done
mv __out-all2 ./$outdir-$CC/out-$MM.pdb
grep $outdir ./$outdir-$CC/out-$MM.pdb > ./$outdir-$CC/out2-$MM.pdb
done

awk '!seen[$0]++' __out-all > ./$outdir-$CC/$outdir-$CC-all.pdb
grep $outdir ./$outdir-$CC/$outdir-$CC-all.pdb > ./$outdir-$CC/$outdir-$CC-all2.pdb
rm __*
