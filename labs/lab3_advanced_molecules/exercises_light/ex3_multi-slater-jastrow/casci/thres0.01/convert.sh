#!/bin/bash

convert4qmc h2o.casci.out -ci h2o.casci.out -readInitialGuess 57 -add3BodyJ > c4q.out
mv sample.Gaussian-G2.xml h2o.wfs.xml
mv sample.Gaussian-G2.ptcl.xml h2o.ptcl.xml
