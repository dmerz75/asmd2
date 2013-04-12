#!/bin/bash
# rsync ggate

function rsy {
    rsync -avh $1 $2
 }
function sc {
    scp -r $1 $2
 }
echo $1
CURRENT=`pwd`

D=hbureau3@fgate-fs.chemistry.gatech.edu:/nethome/hbureau3
S=$CURRENT/$1
sc $S $D
