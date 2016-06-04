#!/bin/sh


for i in $1/*.jpg;
do
    convert $i -crop $2x$3+$4+$5 $(echo $i | cut -d . -f 1).png
    echo "convert $i -crop $2x$3+$4+$5 $(echo $i | cut -d . -f 1).png"
done
