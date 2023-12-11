VMC with zero Jastrows should be equal to HF. 

>> qmca -q eV vmc*/*scalar.dat
                            LocalEnergy               Variance           ratio
vmc_-e/vmc  series 0  -97.010090 +/- 0.009582   20.622493 +/- 0.154381   0.2126
vmc_+e/vmc  series 0  -98.297965 +/- 0.009332   19.529750 +/- 0.118107   0.1987
vmc_optical/vmc  series 0  -97.524116 +/- 0.009116   19.909692 +/- 0.135411   0.2042
vmc/vmc  series 0  -97.790689 +/- 0.008861   21.981607 +/- 1.897721   0.2248

With these results 
QP gap      7.437 +/- 0.498 eV
Optical gap 7.253 +/- 0.345 eV

