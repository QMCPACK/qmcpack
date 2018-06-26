# QMCPACK Manual

This directory contains the LaTeX source for the QMCPACK manual. The
PDF manual is not currently created by the build process and must be
created manually.

- A texlive installation of 1/20/2017 or later included the xetex-collectionis suggested to build the pdf manual.
- An up to date texlive installation is suggested to build the html version of the manual.

---
ATTENTION:
If you have a texlive installation predating the 1/20/2017 release or not including the xetex collection
use the build for legacy tex installs.
``` shell
build_pdflatex_manual_legacy.sh
```
The content is the same but various formatting issues exist. The manual can be built with only the files in the manual directory, consider building it somewhere with a current tex distribution.
---

A script, build_manual.sh provides the current "best" version of the manual.


## OS X
For osx a MacTex of a vintage equivalent equivalent to texlive(1/20/2017) release is required. You will also need to symlink the texlive truetype fonts in so xetex can find them:
```
ln -s /usr/local/texlive/2018/texmf-dist/fonts/truetype/ ~/Library/Fonts/texlive-truetype
```

## Linux
Depending on your method of texlive install you many also need to do minor configuration the get the XCharter true type font included with texlive. Consult documentation for your distro.

# QMCPACK HTML Manual
In addition to the suggestion of up-to-date texlive (do not assume anything is broken unless you have it). `pdf2svg` is required, if it is not available to you via a package on clean Centos7 this is a path to fufilling this requirement.

## Recipe for pdf2svg on clean Centos7 VM
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
```
It is assume you have [spack](https://github.com/spack/spack) to assist in dealing with package management. Yum is used to install packages that should (or must _poppler-glib-devel_) be available at the system level.
