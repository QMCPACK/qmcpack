#!/usr/bin/env bash

if [ "$#" -ne 2 ]; then
    echo "Usage: $0 <first_run_prefix> <restart_run_prefix>"
    echo "Example: $0 qmc_batch.s001 qmc_batch.s002"
    exit 1
fi

file1="${1}.scalar.dat"
file2="${2}.scalar.dat"

echo "========================================================================================"
echo "First run : $file1"
cat "$file1"

echo "========================================================================================"
echo "Restart run : $file2"
cat "$file2"

echo "========================================================================================"
echo "comparing the Kinetic LocalECP NonLocalECP ElecElec IonIon up to the 7th digit after the decimal"

sed "/#/d" "$file1" > s001
sed "/#/d" "$file2" > s002

awk -v tol=1e-6 '
function abs(v) { return v < 0 ? -v : v }
{
  half = NF / 2;
  fail = 0;
  for (col=5; col<=9; col++) {
     val1 = $col;
     val2 = $(col + half);
     diff = abs(val1 - val2) / (abs(val1) + 1e-15);
     if (diff > tol) {
       printf "Column %d differs! val1=%e, val2=%e, rel_diff=%e\n", col, val1, val2, diff;
       fail = 1;
     }
  }
  if (fail) exit 1;
}
' <(paste s001 s002)
