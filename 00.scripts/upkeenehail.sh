#!/bin/bash
# rsync steele

function rsy {
    rsync -avh $1 $2
 }
function sc {
    scp -r $1 $2
 }
echo $1
CURRENT=`pwd`

D=hbureau@keeneland.gatech.xsede.org:/lustre/medusa/hbureau
S=$CURRENT/$1
sc $S $D
