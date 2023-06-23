#!/usr/bin/env bash

echo "========================================================================================"
echo "First run : det_qmc_vmcbatch_dmcbatch_tm.s002.scalar.dat"
cat det_qmc_vmcbatch_dmcbatch_tm.s002.scalar.dat

echo "========================================================================================"
echo "Restart run : det_qmc_vmcbatch_dmcbatch_tm.s003.scalar.dat"
cat det_qmc_vmcbatch_dmcbatch_tm.s003.scalar.dat

echo "========================================================================================"
echo "comparing the Kinetic ElecElec IonIon LocalECP up to the 7th digit after the decimal"

sed "/#/d" det_qmc_vmcbatch_dmcbatch_tm.s002.scalar.dat | awk '{printf("%.7e %.7e %.7e %.7e\n",$5,$6,$7,$8)}' > s002
sed "/#/d" det_qmc_vmcbatch_dmcbatch_tm.s003.scalar.dat | awk '{printf("%.7e %.7e %.7e %.7e\n",$5,$6,$7,$8)}' > s003
diff s002 s003
