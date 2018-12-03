#!/bin/sh

# Legacy script to build PDF version of manual using pdflatex
# For use when xetex is not installed
# Using build_manual.sh with a modern TeX Live installation is the preferred route

if [ $# -eq 1 ]; then
  QMCPACK_VER="$1"
  cp -p version.tex version.save.tex
  sed -i "s/Development Version/$QMCPACK_VER/" version.tex
fi

pdflatex qmcpack_manual.tex
bibtex qmcpack_manual1
bibtex qmcpack_manual2
pdflatex qmcpack_manual.tex
pdflatex qmcpack_manual.tex

if [ ! -z "$QMCPACK_VER" ]; then
  mv version.save.tex version.tex
  echo Added QMCPACK version $QMCPACK_VER
fi
