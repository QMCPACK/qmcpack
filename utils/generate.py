#!/bin/env python
#
# This is a script to generate a series of input files from a template
# and a table of variable values.  The template file can have any number
# of variables of the for ${varname} for replacement.  The data file 
# should then have the following structure:
# filename     V    id    lattice
# file40.xml   40   0      5.31414
# file45.xml   45   1      5.52312
# ...
# This script will then generate file40.xml, file50.xml, ..., replacing
# all instances of ${V}, ${id}, and ${lattice} in the files with the
# value specified for each file.

import sys

if (len(sys.argv) < 3) :
    print "Usage:  generate.py template.txt data.dat"

tempfile = open(sys.argv[1], 'r')
template = tempfile.read()
tempfile.close()

data     = open(sys.argv[2], 'r')
varlist = data.readline().split()[1::]

for line in data:
    fname = line.split()[0]
    args  = line.split()[1::]
    vardict = {}
    for i in range(0,len(varlist)):
        vardict['${'+varlist[i]+'}'] = args[i]
    
    tempcopy = template.replace('a', 'a')
    for key in vardict.keys():
        tempcopy = tempcopy.replace(key, vardict[key])
    
    newfile = open(fname, 'w')
    newfile.write(tempcopy)
    newfile.close()
