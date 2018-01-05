# QMCPACK Manual

This directory contains the LaTeX source for the QMCPACK manual. The
PDF manual is not currently created by the build process and must be
created manually. Note that due to the use of the bibtopic package to
support two bibliographies, bibtex must be invoked twice:

```
pdflatex qmcpack_manual.tex
bibtex qmcpack_manual1
bibtex qmcpack_manual2
pdflatex qmcpack_manual.tex
pdflatex qmcpack_manual.tex
```

