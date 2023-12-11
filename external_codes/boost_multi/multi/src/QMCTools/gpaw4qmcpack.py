#!/usr/bin/env python3
######################################################################################
## This file is distributed under the University of Illinois/NCSA Open Source License.
## See LICENSE file in top directory for details.
##
## Copyright (c) 2021 QMCPACK developers.
##
## File developed by: Juha Tiihonen, tiihonen@iki.fi, University of Jyvaskyla
##
## File created by: Juha Tiihonen, tiihonen@iki.fi, University of Jyvaskyla 
#######################################################################################

# This file converts GPAW orbitals for QMCPACK using ESHDF Python classes

import argparse

from Eshdf_gpaw import EshdfFilePwGpaw

parser = argparse.ArgumentParser(description='GPAW4QMCPACK')
parser.add_argument('infile',help='input .gpw restart file')
parser.add_argument('outfile',default='eshdf.h5',help='output .h5 orbital file')
parser.add_argument('-d','--density',action='store_true',help='whether or not to convert density')

if __name__=='__main__':
    args = vars(parser.parse_args())
    eshdf = EshdfFilePwGpaw(**args)
    eshdf.write()
    print('Converted GPAW orbitals for QMCPACK!')
#end if
