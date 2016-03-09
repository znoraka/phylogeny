#!/bin/sh

cd /usr/include/ImageMagick-6/Magick++

for var in *
do
    echo $var && cat $var | grep ImageStatistics
done
