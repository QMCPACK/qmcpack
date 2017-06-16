#!/bin/bash

echo "========================================================================================"
echo "First run : qmc_short.s001.scalar.dat"
cat qmc_short.s001.scalar.dat

echo "========================================================================================"
echo "Restart run : qmc_short.s002.scalar.dat"
cat qmc_short.s002.scalar.dat

echo "========================================================================================"
echo "comparing the Kinetic ElecElec IonIon LocalECP up to the 7th digit after the decimal"

sed "/#/d" qmc_short.s001.scalar.dat | awk '{printf("%.7e %.7e %.7e %.7e\n",$5,$6,$7,$8)}' > s001
sed "/#/d" qmc_short.s002.scalar.dat | awk '{printf("%.7e %.7e %.7e %.7e\n",$5,$6,$7,$8)}' > s002
diff s001 s002
