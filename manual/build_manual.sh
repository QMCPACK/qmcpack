#!/bin/sh
echo "----♦----♦----♦----♦----♦----♦----♦----♦----♦----♦----♦----♦----"
echo "Building unicode/xelatex (XƎTEX) version of QMCPACK manual."
echo "----♦----♦----♦----♦----♦----♦----♦----♦----♦----♦----♦----♦----"
echo "If you don't see diamonds and a backwards E in the banner"
echo "your terminal and/or editor cannot display unicode characters"
echo "If additionally you don't see a logical representation of the"
echo "unicode code points above in your editor You should show caution"
echo "editing such characters in the .tex files."
echo
if [ -z "$(command -v xelatex)" ]; then
  echo "Error: xelatex is required for build_manual.sh"
  echo
  echo "If typegraphic quality is important to you"
  echo "please install xelatex, preferably as part of"
  echo "texlive(1/20/2017) or later."
  echo "Otherwise if you have a working pdflatex"
  echo "run 'build_pdflatex_manual_legacy.sh'"
  exit 1
fi

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
xelatex qmcpack_manual.tex

if [ ! -z "$QMCPACK_VER" ]; then
  mv version.save.tex version.tex
  echo Added QMCPACK version $QMCPACK_VER
fi
exit 0
