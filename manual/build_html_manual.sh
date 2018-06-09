#!/bin/sh

if [ $# -eq 1 ]; then
  QMCPACK_VER="$1"
  cp -p version.tex version.save.tex
  sed -i "s/Development Version/$QMCPACK_VER/" version.tex
fi
htxelatex qmcpack_manual.tex "xhtml,hmtl5,2,svg,charset=utf8" "-svg -utf8 -v2" "-cvalidate"
bibtex qmcpack_manual1
bibtex qmcpack_manual2
htxelatex qmcpack_manual.tex "xhtml,hmtl5,2,svg,charset=utf8" "-svg -utf8 -v2" "-cvalidate"
htxelatex qmcpack_manual.tex "xhtml,hmtl5,2,svg,charset=utf8" "-svg -utf8 -v2" "-cvalidate"

if [ ! -z "$QMCPACK_VER" ]; then
  mv version.save.tex version.tex
  echo Added QMCPACK version $QMCPACK_VER
fi
