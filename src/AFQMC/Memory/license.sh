
string="///////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source
// License.  See LICENSE file in top directory for details.
//
// Copyright (c) 2020 QMCPACK developers.
//
// File developed by: Fionn Malone, malone14@llnl.gov, LLNL
//
// File created by: Fionn Malone, malone14@llnl.gov, LLNL
////////////////////////////////////////////////////////////////////////////////
"
files=$(find . -iregex '.*\.\(hpp\|h\|cpp\)$')
#files=(blas_hip_catch_all.hpp)
for f in ${files[@]}; do
    echo $f
    sed -i 1,14d $f
    echo "$string" | cat - $f > temp && mv temp $f
done
