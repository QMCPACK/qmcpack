#!/bin/bash

# default threshold is 0.01
convert4qmc h2o.cisd.out -ci h2o.cisd.out -readInitialGuess 57 -add3BodyJ
mv sample.Gaussian-G2.xml h2o.wfs.xml
mv sample.Gaussian-G2.ptcl.xml h2o.ptcl.xml
