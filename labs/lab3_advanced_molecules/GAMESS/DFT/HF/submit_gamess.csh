#!/bin/csh

# location of GAMESS 
set gms = "/soft/applications/qmcpack/DFT_Binaries/gmsbgq"

# do not change
set nodes = 1
set time = 20
set Project = QMC_2014_training
set mode = c2

# name of input file (without .inp)
set name = H2O.HF

$gms ${name} ${nodes} ${time} ${mode} ${Project} ${name}
