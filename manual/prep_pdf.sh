#!/bin/bash

for file in $1; do
    filename=${file%.*}
    pngfile=${filename}.png
    stat $pngfile
    #if the whatever doesn't exist
    if [ $? -ne 0 ]
    then
       echo "converting ${filename}.pdf to ${filename}.png"
       convert -density 300x300 "$filename.pdf" "$filename.png"
    fi
    svgfile=${filename}.svg
    stat $svgfile
    if [ $? -ne 0 ]
    then
       echo "converting ${filename}.pdf to ${filename}.svg"
       pdfcrop --margins 10 --clip --resolution 144 "$filename.pdf" "$filename.pdf"
       pdf2svg "$filename.pdf" "$filename.svg"
       ebb -x ${filename}.pdf
    fi
done
