#! /usr/bin/env bash

if [ ! -e "$MXDIR/scripts/paper" ]; then
    echo "Please invoke using correct MXDIR as MXDIR=$MXDIR $MXDIR/fixpaths.sh"
    exit
fi

cd $MXDIR/scripts/paper

for f in $(\ls 0*.sh); do sed -i 's;^MMFITDIR=~/.*$;MMFITDIR=${MXDIR}/src;' $f; done

sed -i 's;^ORIGMDIR=\$W.*;ORIGMDIR=${MXDIR}/data/numtests;' 007gendat.sh
sed -i 's;^ORIGCOVDIR=\$W.*;ORIGCOVDIR=${MXDIR}/data/numtests;' 007gendat.sh
sed -i 's;^ORIGMDIR=\$W.*;ORIGMDIR=${MXDIR}/data/numtests;' 043gendat.sh
sed -i 's;^ORIGCOVDIR=\$W.*;ORIGCOVDIR=${MXDIR}/data/numtests;' 043gendat.sh

for f in $(\ls 0*.sh); do sed -i 's;/vol/3/mmfit2022;${MXDIR};' $f; done



for f in $(\ls [34]*.sh); do sed -i 's;^MMFITDIR=~/.*$;MMFITDIR=${MXDIR}/src;' $f; done

sed -i 's;^DATADIR=/vol/.*;DATADIR=${MXDIR}/data/cyp2/perchain;' 382group.sh
sed -i 's;^DATADIR=/vol/.*;DATADIR=${MXDIR}/data/abl/perchain;' 411ablgroup.sh
sed -i 's;^DATADIR=/vol/.*;DATADIR=${MXDIR}/data/cam/perchain;' 612camgrpc145nx1.sh

for f in $(\ls [346]*.sh); do sed -i 's;/vol/3/mmfit2022;${MXDIR};' $f; done

