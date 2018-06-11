#!/bin/bash
set -x
for file in $1; do
    filename=${file%.*}
    echo "converting ${filename}.pdf to ${filename}.svg"
    pdfcrop --margins 10 --clip --resolution 144 "$filename.pdf" "$filename.pdf"
    pdf2svg "$filename.pdf" "$filename.svg"
    echo "converting ${filename}.pdf to ${filename}.png"
    convert -density 300x300 "$filename.pdf" "$filename.png"
done
