#! /usr/bin/env bash

DATADIR="$1"

destdir=./
classfile=clustmap.txt

subtypes="$(for f in $(ls $DATADIR/*.pdb); do basename $f | sed 's/-.*//'; done | sort | uniq)"

for st in ${subtypes}; do
    cp /dev/null ${st}.pdb
done    
  ## redefine subtypes(enzymes) according to files on $destdir
subtypes="$(for f in $(ls $destdir/*.pdb); do basename $f .pdb | sed 's/-.*//'; done | sort | uniq)"

echo $subtypes

cp /dev/null pdbmfiles.txt
for st in ${subtypes}; do
    mi=0
    cp /dev/null ${st}.pdb
    for f in $(ls -d $DATADIR/${st}-*.pdb); do
	((mi++))
	printf "MODEL %8d\n" ${mi} >> ${st}.pdb
	cat $f >> ${st}.pdb
	echo "ENDMDL" >> ${st}.pdb
    done
    echo "END" >> ${st}.pdb
done

cp /dev/null pooled.pdb
mi=0
k=0
echo pdbid cid > $classfile
for st in ${subtypes}; do
    for f in $(ls -d $DATADIR/${st}-*.pdb); do
	((mi++))
	printf "MODEL %8d\n" ${mi} >> pooled.pdb
	cat $f >> pooled.pdb
	echo "ENDMDL" >> pooled.pdb
        echo \"$(basename $f .pdb)\" $k >> $classfile
    done
    ((k++))
done
echo "END" >> pooled.pdb
