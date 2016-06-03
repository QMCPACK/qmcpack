#!/bin/bash

convert4qmc h2o.cas_rerun.out -ci h2o.cas_rerun.out -threshold 0.0075 -add3BodyJ > c4q.out
mv sample.Gaussian-G2.xml h2o.wfs.xml
mv sample.Gaussian-G2.ptcl.xml h2o.ptcl.xml
