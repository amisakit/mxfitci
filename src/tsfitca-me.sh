#!/usr/bin/env bash

MMFITDIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )

MSFIT=$MMFITDIR/msfit
MARKPDB=$MMFITDIR/markpdb

NENSLIM=99999

## gets MEDOIDs from d1${mols}.bin and write them on ts-ave-${mols}.bin
# $1: natoms
# $2: mols
# $3: medoid k, 1 <= k <= mi
R_write_ens_aves_bin() {
    R --vanilla --slave --args $1 $2 $3 <<EOF
    args <- commandArgs(trailing=T)
    natoms <- as.numeric(args[1])
    mols <- args[2]
    k <- as.numeric(args[3])
    con <- file(paste0('d1',mols,'.bin'),'rb')
    pos <- seek(con, where=(natoms^2 + (3*natoms)*k)*8,origin='start')
    m <- readBin(con, 'double', 3*natoms)
    close(con)
    con <- file(paste0('ts-ave-',mols,'.bin'),'wb')
    writeBin(m, con)
    close(con)
EOF
}

## read V, M, V9, and Rbuf from `dts.bin`,
## then write them on respective files.
# $1: natoms
# $2: nens - numaborted (effective nens)
R_write_MVR_bin() {
    R --vanilla --slave --args $1 $2 <<EOF
    args <- commandArgs(trailing=T)
    natoms <- as.numeric(args[1])
    nens <- as.numeric(args[2])
    con <- file('dts.bin','rb')
    v <- readBin(con, 'double', natoms^2)
    m <- readBin(con, 'double', 3*natoms)
    seek(con, where=(natoms^2 + 3*natoms + nens*3*natoms)*8, origin='start')
    v9 <- readBin(con, 'double', 9*natoms^2)
    rbuf <- readBin(con, 'double', 9*nens+3*nens)
    close(con)
    con <- file('tsV.bin', 'wb')
    writeBin(v, con)
    close(con)
    con <- file('m.bin','wb')
    writeBin(m, con)
    close(con)
    con <- file('tsV9.bin','wb')
    writeBin(v9, con)
    close(con)
    con <- file('tsrbuf.bin','wb')
    writeBin(rbuf, con)
    close(con)
EOF
}

## mi-weighted average of Sj's
## variable mi's are supported (not yet verified).
## nens, natoms, and mi's are get from `i1.*.bin`.
## Rbuf are read from `tsrbuf.bin`, Sj and S9j are read from `d1.*.bin`.
# arguments: none
R_weighted_ave_S() {
    R --vanilla --slave <<EOF
    a <- lapply(list.files(pattern='i1.*.bin'), function(x) {con <- file(x,'rb'); nm <- readBin(con, 'int', 2); close(con); nm})
    nens = length(a)
    natoms <- a[[1]][1]
    mis <- array(as.numeric(unlist(a)), c(length(a[[1]]),length(a)))[2,]

    con <- file('tsrbuf.bin', 'rb')
    rbuf <- array(readBin(con, 'double', 9*nens), c(3,3,nens))
    close(con)

    dfiles <- list.files(pattern='d1.*.bin')

    a <- lapply(dfiles, function(x) {con <- file(x,'rb'); v <- readBin(con, 'double', natoms^2); close(con); v})
    con <- file('tsS.bin', 'wb')
    writeBin(apply(do.call('cbind', a), 1, function(x) weighted.mean(x, mis)), con)
    close(con)

    a <- lapply(1:nens, function(k) {con <- file(dfiles[k],'rb'); seek(con, where=(natoms^2 + 3*natoms + mis[k]*3*natoms)*8, origin='start'); s9 <- readBin(con, 'double', 9*natoms^2); close(con); s9})
    a <- do.call('cbind', a)
    b <- aperm(array(a,c(3,natoms,3,natoms,nens)), c(1,3,2,4,5))
    br <- lapply(1:nens, function(k) {R <- rbuf[,,k]; apply(b[,,,,k], c(3,4), function(x) { R %*% x %*% t(R) })})
    b <- aperm(array(do.call('cbind',br), c(3,3,natoms,natoms,nens)), c(1,3,2,4,5))
    a <- apply(b, c(1,2,3,4), function(x) weighted.mean(x, mis))

    con <- file('tsS9.bin', 'wb')
    writeBin(as.numeric(a), con)
    close(con)

    # Siso by averaging Saniso (tsS9)
    b <- aperm(array(a,c(3,natoms,3,natoms)), c(1,3,2,4))
    a <- apply(b, c(3,4), function(x) {mean(diag(x))})
    con <- file('tsSave.bin', 'wb')
    writeBin(as.numeric(a), con)
    close(con)
EOF
}

# do.call('rbind', a) == t(array(as.numeric(unlist(a)),c(nrow,ncol)))
# do.call('cbind', a) == array(as.numeric(unlist(a)),c(nrow,ncol))


estimate(){
    # ensemble averages
    aves=()
    numaborted=0
    for k in `seq 0 $((nens-1))`; do
        mols=a$(printf "%05d" $k)                # NENSLIM
        pdbm=${pdbmfiles[$k]}
	$MSFIT $WEIGHT -v d -e 1 1 1 1 500 -o c _tmp1 -o s ts-cov-$mols.dat -b o 1$mols.bin $pdbm 2>&1 | tee msfit1-$mols.log
	if [ -f d1$mols.bin ]; then
	    sed -e '/====== mean/,$d' -e 's/  *@.*//' _tmp1 > ts-crd-$mols.pdb
	    k=$(cat msfit1-$mols.log \
		|  awk 'BEGIN {prt=0} /^$/ { prt=0 }  {if (prt==1) {print $5}} /^   id/ { prt=1 }' | R --slave -e 'cat(sort(scan(file="stdin",quiet=T),index.return=T,decreasing=F)$ix[1],"\n")')
	    R_write_ens_aves_bin $natoms $mols $k
	    cat _tmp1 \
		| awk -v modelid=$k 'BEGIN {prt=0; atom=0} /^MODEL/ { ++prt; atom=1 } {if (prt==modelid && atom==1) {print $0}} /^ENDMDL/ { atom=0 }' | grep '^ATOM' > _tmp2
	    echo "TER" >> _tmp2
	    $MARKPDB "CA" _tmp2 > ts-ave-$mols.pdbm
	    aves+=( ts-ave-$mols.pdbm )
	else
	    (( numaborted++ ))
	fi
	rm -f _tmp1 _tmp2
    done

    # grand average
    $MSFIT $WEIGHT -v d -e 1 1 1 1 500 -o c _tmp3 -o s tsV.dat -b i bin -b o ts.bin ${aves[@]} 2>&1 | tee msfit2.log
    R_write_MVR_bin $natoms $(( nens - numaborted ))
    sed 's/  *@.*//' _tmp3 > ts-aves.pdb		# strip landmarks
    $MARKPDB "CA" ts-aves.pdb > _tmp4			# reassign landmarks
    rm -f ts-aves.pdb _tmp3

    i=$(( (natoms+1)*(nens+1-numaborted) ))
    if [ `wc -l _tmp4 | awk '{print $1}'` -ne $i ]; then
	echo "(natoms = $natoms + 1 for TER) x ($nens + 1 - $numaborted) = $i lines are assumed."
	exit -1
    fi

    # strip landmarks on the ensemble averages keeping those on the grand ave
    i=$(( i - natoms - 1 ))
    cat _tmp4 | sed "1,${i}s/  *@.*//" > ts-aves.pdbm

    # pooled sample covariance as an estimate for Sigma
    R_weighted_ave_S
    echo $numaborted > numaborted.txt
    rm -f _tmp4
    return
}


R_append_mis() {
    R --vanilla --slave <<EOF
con <- file('its.bin','rb')
natoms <- readBin(con,'int')
nens <- readBin(con,'int')
close(con)
mis <- unlist(lapply(list.files(pattern='i1a0.*.bin'), function(x) {con <- file(x,'rb'); readBin(con,'int'); mi <- readBin(con,'int'); close(con); mi}))
if (length(mis) != nens) { cat("inconsistent mi's and nens\n"); q(status=-1) }
con <- file('its.bin','ab')
writeBin(mis, con)
close(con)
q(status=0)
EOF
    if [ $? -ne 0 ]; then
        echo Error
        exit -1
    fi
}


# tsfitca.sh [-0] $1
# $1: a file containing the list of pdbm files (full paths)
# -0: OLS fit

usage() { echo "$(basename $0) [-0] pdbmlist.txt"; exit; }

WEIGHT=

while getopts "0" opt; do
    case $opt in
        0) WEIGHT=-h
           ;;
        *) usage
           ;;
    esac
done
shift $((OPTIND-1))

pdbmfiles=(`cat $1`)
nens=${#pdbmfiles[@]}
if [ $nens -gt $NENSLIM ]; then
  echo "r too large"
  exit
fi
mi1=`grep MODEL ${pdbmfiles[1]} | wc -l`
natoms=$((`grep @ ${pdbmfiles[1]} | wc -l` / mi1))

 echo natoms=$natoms
 echo nens=$nens
 echo mi=$mi1
 echo pdbmfiles=${pdbmfiles[1]},...

estimate

#rm ts-ave-*.{bin,pdbm} ts-cov-*.dat ts-crd-*.pdb i1a*.bin d1a*.bin

R_append_mis
