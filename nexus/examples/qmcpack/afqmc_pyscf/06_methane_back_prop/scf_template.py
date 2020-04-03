#! /usr/bin/env python3

from math import cos, sin, pi, acos
import numpy
from pyscf import gto, scf

$system

mf = scf.RHF(mol)
mf.chkfile = $chkfile
ehf = mf.scf()
