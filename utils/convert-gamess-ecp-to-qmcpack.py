#!/usr/bin/env python3

import sys,os,getopt
sys.path.append(os.path.join(os.path.dirname(__file__), '../nexus/lib'))

executable_name = sys.argv[0]
inputfile = None
try:
    opts, args = getopt.getopt(sys.argv[1:],"hi:",["ifile="])
except getopt.GetoptError:
    print('Usage:\n{} -i <gamess-input-pp>'.format(executable_name))
    sys.exit(2)
for opt, arg in opts:
    if opt == '-h':
        print('Usage:\n{} -i <gamess-input-pp>'.format(executable_name))
        sys.exit()
    elif opt in ("-i", "--ifile"):
        inputfile = arg
#end try

if inputfile is None:
    print('Usage:\n{} -i <gamess-input-pp>'.format(executable_name))
    sys.exit(2)
#end if

from pseudopotential import GaussianPP
# Read semi-local potential
pp = GaussianPP(inputfile,format='gamess')
# write to qmcpack format
with open('{}.xml'.format(inputfile),'w') as f:
    f.write(pp.write_text(format='qmcpack'))
#end with
