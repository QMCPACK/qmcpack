# QMCPACK Manual

This directory contains the LaTeX source for the QMCPACK manual. The
PDF manual is not currently created by the build process and must be
created manually. Note that due to the use of the bibtopic package to
support two bibliographies, bibtex must be invoked twice:

```
xelatex qmcpack_manual.tex
bibtex qmcpack_manual1
bibtex qmcpack_manual2
xelatex qmcpack_manual.tex
xelatex qmcpack_manual.tex
```

A script, build_manual.sh, is provided.

