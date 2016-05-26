#!/bin/sh

COUNTER=0

for f in $1/*.pgm
do
    echo "file = " $f
    mkdir $2/$COUNTER
    
    ../../tree-generator/bin/tree-generator "/media/ramdisk/data" $f 10 0 ./$COUNTER
    ../../tree-extractor/c++/bin/tree-extractor "/media/ramdisk/data/" ./$COUNTER
    ../../tree-comparator/bin/tree-comparator ./$COUNTER/truth.txt ./$COUNTER/computed.txt > ./$COUNTER/results.txt
    
    COUNTER=$((COUNTER + 1))
done
