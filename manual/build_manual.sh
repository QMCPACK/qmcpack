#!/bin/sh

if [ $# -eq 1 ]; then
  QMCPACK_VER="$1"
  cp -p version.tex version.save.tex
  sed -i "s/Development Version/$QMCPACK_VER/" version.tex
fi

xelatex qmcpack_manual.tex
bibtex qmcpack_manual1
bibtex qmcpack_manual2
xelatex qmcpack_manual.tex
xelatex qmcpack_manual.tex

if [ ! -z "$QMCPACK_VER" ]; then
  mv version.save.tex version.tex
  echo Added QMCPACK version $QMCPACK_VER
fi
