# QMCPACK Manual

This directory contains the LaTeX source for the QMCPACK manual. The
PDF manual is not currently created by the build process and must be
created manually.

A texlive installation of 2017 or later included xetex packages is suggested to build the pdf manual.

Note that due to the use of the bibtopic package to
support two bibliographies, bibtex must be invoked twice:

```
xelatex qmcpack_manual.tex
bibtex qmcpack_manual1
bibtex qmcpack_manual2
xelatex qmcpack_manual.tex
xelatex qmcpack_manual.tex
```

A script, build_manual.sh, is provided.


## OS X
For osx a MacTex of a vintage equivalent to texlive is required.
You will also need to symlink the texlive truetype fonts in so xetex can find them:
```
ln -s /usr/local/texlive/2018/texmf-dist/fonts/truetype/ ~/Library/Fonts/texlive-truetype
```

## Old or limited tex installations
```
pdflatex qmcpack_manual.tex
bibtex qmcpack_manual1
bibtex qmcpack_manual2
pdflatex qmcpack_manual.tex
pdflatex qmcpack_manual.tex
```

A script, build_pdflatex_manual.sh, is provided.


