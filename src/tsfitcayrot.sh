#!/usr/bin/env bash

MMFITDIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )

MSFIT=$MMFITDIR/msfit
MARKPDB=$MMFITDIR/markpdb

NENSLIM=99999

#
# rotate Mi and Yij according to the second stage rotation Ri,
# in the current directory
# arguements: none
R_rotate_Y() {
    R --vanilla --slave <<EOF

    a <- lapply(list.files(pattern='i1.*.bin'), function(x) {con <- file(x,'rb'); nm <- readBin(con, 'int', 2); close(con); nm})
    nens = length(a)
    natoms <- a[[1]][1]
    mis <- array(as.numeric(unlist(a)), c(length(a[[1]]),length(a)))[2,]

    con <- file('tsrbuf.bin', 'rb')
    rbuf <- array(readBin(con, 'double', 9*nens), c(3,3,nens))
    tbuf <- array(readBin(con, 'double', 3*nens), c(3,nens))
    close(con)

    dfiles <- list.files(pattern='d1.*.bin')

    out <- file('yrot.bin', 'wb')
    mylist <- lapply(1:nens, function(k) {
        con <- file(dfiles[k],'rb');
        seek(con, where=(natoms^2)*8, origin='start');
        my <- array(readBin(con, 'double', 3*natoms*(1+mis[k])), c(3,natoms*(1+mis[k])));
        close(con);
        a <- rbuf[,,k] %*% (my - tbuf[,k]);
        writeBin(as.numeric(a), out);
    })
    close(out)    
EOF
}

R_confirm_rotation() {
    R --vanilla --slave <<EOF

    a <- lapply(list.files(pattern='i1.*.bin'), function(x) {con <- file(x,'rb'); nm <- readBin(con, 'int', 2); close(con); nm})
    nens = length(a)
    natoms <- a[[1]][1]
    mis <- array(as.numeric(unlist(a)), c(length(a[[1]]),length(a)))[2,]

    con2 <- file('dts.bin', 'rb')
    pos <- seek(con2, (natoms^2 + 3*natoms) * 8, origin='start')
    con1 <- file('yrot.bin', 'rb')
    dum <- lapply(1:nens, function(k) {
        m1 <- readBin(con1, 'double', 3*natoms);
        pos <- seek(con1, (3*natoms * mis[k]) * 8, origin="current");
        m2 <- readBin(con2, 'double', 3*natoms);
        cat('err', sum((m1-m2)^2), '\n');
    })
    close(con1)
    close(con2)
EOF
}


R_rotate_Y

R_confirm_rotation
