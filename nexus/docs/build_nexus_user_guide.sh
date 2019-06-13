#!/bin/sh
echo "----♦----♦----♦----♦----♦----♦----♦----♦----♦----♦----♦----♦----"
echo "Building unicode/xelatex (XƎTEX) version of NEXUS User Guide."
echo "----♦----♦----♦----♦----♦----♦----♦----♦----♦----♦----♦----♦----"
echo "If you don't see diamonds and a backwards E in the banner"
echo "your terminal and/or editor cannot display unicode characters"
echo "If additionally you don't see a logical representation of the"
echo "unicode code points above in your editor You should show caution"
echo "editing such characters in the .tex files."
echo
if [ -z "$(command -v xelatex)" ]; then
  echo "Error: xelatex is required for build_nexus_user_guide.sh"
  exit 1
fi

# QMCPACK build_manual.sh has examples on versioning .tex and .pdf
# NEXUS version is currently fixed

xelatex -shell-escape nexus_user_guide.tex
xelatex -shell-escape nexus_user_guide.tex
xelatex -shell-escape nexus_user_guide.tex

exit 0
