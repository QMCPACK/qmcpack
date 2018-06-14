#!/bin/sh
# Builds the html manual
# texlive install suggested

if [ $# -eq 1 ]; then
  QMCPACK_VER="$1"
  cp -p version.tex version.save.tex
  sed -i "s/Development Version/$QMCPACK_VER/" version.tex
fi

cp -r figures html_site/
./prep_pdf.sh './html_site/figures/*.pdf'
make4ht -u -x -e qmcpack_manual.mk4 -c qmcpack_manual.cfg qmcpack_manual.tex "xhtml,2,html5,graphics-144" "" " -cvalidate" -d "html_site"

if [ ! -z "$QMCPACK_VER" ]; then
  mv version.save.tex version.tex
  echo Added QMCPACK version $QMCPACK_VER
fi
