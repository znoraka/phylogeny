#!/bin/sh

FILES=$1/*.jpg
rm -rf $1/dcts
mkdir $1/dcts
for f in $FILES
do
    # convert -quality 100 $f out.jpg
    x=$(basename $f)
    y=${x%%.*}
    echo $f
    echo $y
    ./dctextractor $f > $1/dcts/$y.txt
    # rm out.jpg
done
