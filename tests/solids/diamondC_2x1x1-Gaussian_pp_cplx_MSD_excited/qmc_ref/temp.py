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

