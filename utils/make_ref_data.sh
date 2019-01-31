#/bin/bash

# This file is distributed under the University of Illinois/NCSA Open Source License.
# See LICENSE file in top directory for details.
#
# Copyright (c) 2019 QMCPACK developers.
#
# File Developed By: Peter Doak, doakpw@ornl.gov, Oak Ridge National Laboratory
#
# File Created By: Peter Doak, doakpw@ornl.gov, Oak Ridge National Laboratory

#this script will generate values to past into the CMakeLists.txt
#for a test against reference with the error_range/sigma
#that check_scalars.py expects

#watch out for ionion which usually has no errorbar and you must pic something by hand
#usually 0.0001

if [ -z $QMCPACK_ROOT ]; then
    QMCPACK_ROOT=/scratch/epd/qmcps/qmcpack
fi

QMCA="${QMCPACK_ROOT}/nexus/bin/qmca"
if [ ! -f $QMCA ]; then
    echo "You must define QMCPACK_ROOT to the directy containing nexus/bin and tests"
fi

if [ "$#" -lt 1 ]; then
    echo "usage: make_ref_data.sh SCALAR_TAG run_file ref_file"
    exit
fi

SCALAR_TAG=$1
if [ "$#" -ge 2 ]; then
    run_file=$2
else
    run_file=qmc_short.s000.scalar.dat
fi
if [ "$#" -eq 3 ]; then
    ref_file=$3
else
    ref_file=qmc-ref/qmc_ref_long.s000.scalar.dat
fi

TEST_PARENT=$(dirname $PWD)
echo `pwd`
echo ${run_file}
echo ${ref_file}
${QMCA} --fp=24.12f -e 3 $run_file 2>/dev/null | gawk -e '/\s\s[a-zA-Z]+\s+=/ {print $1 " " $3 " " $5}' > run_vals
${QMCA} --fp=24.12f -e 3 $ref_file 2>/dev/null | gawk -e '/\s\s[a-zA-Z]+\s+=/ {print $1 " " $3 " " $5}' > ref_vals

PYTHONFILE=temp.py
(
cat <<'EOF'
import numpy as np
import math
import sys

SCALAR_TAG=sys.argv[1]
why_why_why = { "Variance" :"variance",
		"Kinetic": "kinetic",
		"LocalPotential" : "potential",
		"ElecElec" : "eeenergy",  
		"LocalECP" : "localecp",
		"NonLocalECP" : "nonlocalecp",
		"IonIon" : "ionion",
                "LocalEnergy" : "totenergy" }
f = open('run_vals','r')
slines = [ line.split(' ') for line in f.readlines() ]
f = open('ref_vals','r')
rlines = [ line.split(' ') for line in f.readlines() ]

ref_data = []
for ref, short in zip(rlines, slines):
    combined = math.sqrt(np.float64(ref[2])**2 + np.float64(short[2])**2)
    if ref[0] in why_why_why:
        ref_data.append([why_why_why[ref[0]], np.float64(ref[1]), np.float64(combined)])
print(ref_data)
for rd in ref_data:
    print ("LIST(APPEND {0}_SCALARS \"{1}\" \"{2:.12f} {3:.12f}\")".format(SCALAR_TAG, rd[0], rd[1], rd[2]))

EOF
) > $PYTHONFILE

python $PYTHONFILE $SCALAR_TAG
