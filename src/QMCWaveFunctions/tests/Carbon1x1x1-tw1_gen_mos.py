#! /usr/bin/env python3
import numpy
import h5py
from pyscf.pbc import gto, scf, dft, df
from pyscf import __version__ 
import datetime

# Author: Chandler Bennett
# This file is modified from tests/solids/diamondC_1x1x1-Gaussian_pp_Tw_cplx/dft-inputs/Carbon1x1x1-tw1.py
# to print (in C++ syntax) the real and imaginary parts of all molecule orbitals along a path that spans 
# regions both inside and outside the primitive cell. This output was copied into test_pyscf_complex_MO.cpp,
# a unit test that verifies the real/imag values of QMCPACK's molecular orbitals match along the same path.

cell = gto.Cell()
cell.a             = '''
         3.37316115       3.37316115       0.00000000
         0.00000000       3.37316115       3.37316115
         3.37316115       0.00000000       3.37316115'''
cell.atom = '''  
   C        0.00000000       0.00000000       0.00000000
   C        1.686580575      1.686580575      1.686580575 
            ''' 
cell.basis         = 'bfd-vdz'
cell.ecp           = 'bfd'
cell.unit          = 'B'
cell.drop_exponent = 0.1
cell.verbose       = 5
cell.charge        = 0
cell.spin          = 0
cell.build()


sp_twist=[0.333,0.333,0.333]
twist = numpy.asarray([0.333,0.333,0.333]) / 1.0
kmesh=[1,1,1]
kpts = cell.make_kpts((1,1,1), with_gamma_point=False,  wrap_around=True, scaled_center=twist)


mf = scf.KRHF(cell,kpts)
mf.exxdiv = 'ewald'
mf.max_cycle = 200

e_scf=mf.kernel()

import numpy as np
coords=np.zeros((1,3))

print('    // BEGIN generated C++ input from %s (pyscf version %s) on %s\n'%(__file__,__version__,datetime.datetime.now()))
for i in range(7):

    coords[0][0]=-10+3.333333*i
    coords[0][1]=-10+3.333333*i
    coords[0][2]=-10+3.333333*i

    print('    //Move electron 0 to position %s a.u.:'%(str(coords)))
    print('    elec.R[0] = { %s, %s, %s };'%(str(coords[0][0]),str(coords[0][1]),str(coords[0][2])))
    print('    elec.update();')
    print('    sposet->evaluate(elec, 0, values);\n')

    print('    // Position %s a.u.'%(str(coords)))
    ao = cell.pbc_eval_gto('GTOval', coords, kpt=kpts[0])

    for spo_idx in range(len(mf.mo_coeff[0])):

        print('    // Verifying values of SPO %s'%(str(spo_idx)))
        print('    CHECK(std::real(values[%d]) == Approx(%s));'%(spo_idx,str(numpy.dot(numpy.array(mf.mo_coeff)[0].T,ao[0].T)[spo_idx].real)))
        print('    CHECK(std::imag(values[%d]) == Approx(%s));\n'%(spo_idx,str(numpy.dot(numpy.array(mf.mo_coeff)[0].T,ao[0].T)[spo_idx].imag)))

print('    // END generated C++ input from %s (pyscf version %s) on %s'%(__file__,__version__,datetime.datetime.now()))

