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

D=dmerz3@tg-steele.purdue.teragrid.org:/scratch/scratch96/d/dmerz3/valiant/steele/00.interim/
S=$CURRENT/$1
sc $S $D
/media/12e650a0-b7be-4539-a9c4-785021d86825/home/dale/Documents/md/asmd2/00.scripts
