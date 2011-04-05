#!/usr/bin/env bash

rm *.jpg
TAG='70GHz_9x9'
for f in rcond*.png
do
    echo $f
    convert $f $f.jpg
done

montage S*.png -geometry 400x+0+0 Smaps_$TAG.jpg

montage DPC* IQUSS* WMAP* -geometry 400x+0+0 -tile 2x3 QU_$TAG.jpg
