#!/bin/sh

if [ $# -eq 1 ]; then
  QMCPACK_VER="$1"
  cp -p version.tex version.save.tex
  sed -i "s/development/$QMCPACK_VER/" version.tex
fi

pdflatex qmcpack_manual.tex
bibtex qmcpack_manual1
bibtex qmcpack_manual2
pdflatex qmcpack_manual.tex
pdflatex qmcpack_manual.tex

if [ ! -z "$QMCPACK_VER" ]; then
  mv version.save.tex version.tex
  echo added QMCPACK version $QMCPACK_VER
fi
