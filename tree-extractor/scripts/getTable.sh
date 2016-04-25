#!/bin/sh

djpeg -verbose -verbose $1 > /dev/null
echo "Quality =" $(identify -verbose $1 | grep "Quality" | cut -d " " -f4)
