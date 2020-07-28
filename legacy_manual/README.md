# QMCPACK Manual

This directory contains the LaTeX source for the QMCPACK manual. The
PDF manual is not currently created by the build process and must be
created manually. An HTML version can also be created, using
additional tools.

- A script, build_manual.sh provides the current "best" version of the manual.
- A TeX Live installation of 2017-01-20 or later including the xetex-collection is suggested to build the PDF manual.
- An up to date TeX Live installation is suggested to build the HTML version of the manual.

``` shell
build_manual.sh
```

---
# ATTENTION:

If you have a TeX Live installation predating the 2017-01-20 release or not including the xetex collection
use the build for legacy tex installs. The content is the same but various formatting issues exist. Since the manual can be built with only the files in the manual directory, consider building it somewhere with an up to date tex distribution.

``` shell
build_pdflatex_manual_legacy.sh
```
---
# Contributing:
See the chapter "Contributing to the Manual" in contributing.tex
--

## OS X
For OS X, MacTex equivalent to TeX Live(1/20/2017) release or later is suggested. You should also symlink the TeX Live truetype fonts in so xetex can find them:
```
ln -s /usr/local/TeX Live/2018/texmf-dist/fonts/truetype/ ~/Library/Fonts/texlive-truetype
```

# QMCPACK HTML Manual

In addition to the suggestion of up-to-date TeX Live (do not assume
anything is broken unless you have it). `pdf2svg` is required, your system package manager may supply it, homebrew on OSX, and spack on linux are also options.

## Recipes for pdf2svg and other tools

### OSX
Assuming don't have an up to date texlive and for some reason haven't installed homebrew.  If you are still using macports I strongly suggest moving your install tree and giving [homebrew](https://brew.sh) a try.
```shell
/usr/bin/ruby -e "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/master/install)"
brew cask (re)install mactex
brew install pdf2svg
tlmgr update -all
```

### Centos7
It is assumed you have [spack](https://github.com/spack/spack) to assist in dealing with package management. Yum is used to install packages that should (or must _poppler-glib-devel_) be available at the system level.

``` shell
    yum install gcc-c++
    yum install environment-modules
    yum install bzip2
    spack install curl
    spack load curl
    sudo yum install libstdc++-static.x86_64
    spack install gcc@7.3.0 +binutils +piclibs
    spack load gcc@7
    sudo yum install cairo-devel.x86_64
    sudo yum install poppler-glib-devel.x86_64
    spack install pdf2svg%gcc@7.3.0
    tlmgr update make4ht
```
