#!/bin/sh

COUNTER=$(cd $2 && ls -l | grep -v ^l | wc -l)
COUNTER=$((COUNTER - 1))

for i in $(seq 0 $4);
do
    for f in $1/*.pgm
    do
	echo "file = " $f
	mkdir $2/$COUNTER
	
	../../tree-generator/bin/tree-generator "/media/ramdisk/data" $f $3 0 $2$COUNTER
	../../tree-extractor/c++/bin/tree-extractor "/media/ramdisk/data/" $2$COUNTER
	../../tree-comparator/bin/tree-comparator $2$COUNTER/truth.txt $2$COUNTER/computed.txt >> $2$COUNTER/results.txt
	
	COUNTER=$((COUNTER + 1))
    done
done
