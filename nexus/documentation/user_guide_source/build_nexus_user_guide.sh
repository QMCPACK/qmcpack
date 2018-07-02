#!/bin/sh
if [ -z "$(command -v pdflatex)" ]; then
  echo "Error: pdflatex is required for build_nexus_user_guide.sh"
  exit 1
fi

# From QMCPACK build_manual.sh, for eventual versioning
#if [ $# -eq 1 ]; then
#  QMCPACK_VER="$1"
#  cp -p version.tex version.save.tex
#  sed -i "s/Development Version/$QMCPACK_VER/" version.tex
#fi


pdflatex -shell-escape nexus_user_guide.tex
pdflatex -shell-escape nexus_user_guide.tex
pdflatex -shell-escape nexus_user_guide.tex

# From QMCPACK build_manual.sh, for eventual versioning
#if [ ! -z "$QMCPACK_VER" ]; then
#  mv version.save.tex version.tex
#  echo Added QMCPACK version $QMCPACK_VER
#fi
exit 0
