#!/bin/sh
pdflatex qmcpack_manual.tex
bibtex qmcpack_manual1
bibtex qmcpack_manual2
pdflatex qmcpack_manual.tex
pdflatex qmcpack_manual.tex
