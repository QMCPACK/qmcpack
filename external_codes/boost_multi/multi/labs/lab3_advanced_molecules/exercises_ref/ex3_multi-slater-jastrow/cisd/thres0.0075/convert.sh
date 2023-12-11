#!/bin/bash

convert4qmc h2o.cisd.out -ci h2o.cisd.out -readInitialGuess 57 -threshold 0.0075 -add3BodyJ > c4q.out
mv sample.Gaussian-G2.xml h2o.wfs.xml
mv sample.Gaussian-G2.ptcl.xml h2o.ptcl.xml
