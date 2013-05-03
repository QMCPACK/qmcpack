#!/bin/env bash

# shared input files for gr3x3x1 and gr4x4x1
wget http://cms.mcc.uiuc.edu/~jnkim/graphite-benchmark.tgz
tar zxvf graphite-benchmark.tgz

# get hdf5 file for water
cd water
wget http://cms.mcc.uiuc.edu/~jnkim/qmcdb/h2o.pwscf.h5
cd -

cd gr3x3x1
ln -s ../graphite-benchmark/lda.pwscf.h5 .
cd -

cd gr4x4x1
ln -s ../graphite-benchmark/lda.pwscf.h5 .
cd -

