#!/bin/sh
if [ $# -eq 0 ]; then
    echo "No filename supplied"
    exit 1
fi
export MOFDIR=/archive_data/desarchive/OPS/multiepoch/Y3A2_MOF
export LINKDIR=./tilelinks
count=0
if [ ! -d "$LINKDIR" ]; then
    mkdir $LINKDIR
fi
while IFS='' read -r line || [[ -n "$line" ]]; do
    let count=$count+1
    IFS=" " read -ra arr <<< "$line"
    tile=${arr[0]}
    reqnum=${arr[1]}
    attnum=${arr[2]}
    #echo $count $tile
    if [ "$attnum" -lt 10 ]; then 	
        FILE=$MOFDIR/r$reqnum/$tile/p0$attnum/mof/$tile'_r'$reqnum'p0'$attnum'_mof.fits'
    else
        FILE=$MOFDIR/r$reqnum/$tile/p$attnum/mof/$tile'_r'$reqnum'p'$attnum'_mof.fits'
    fi
    if [ -f "$FILE" ]; then
        cd $LINKDIR
        ln -s $FILE
        cd ..
    else
        echo "File $FILE doesn't exist"
    fi
done < "$1"
