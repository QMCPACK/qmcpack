#!/usr/bin/env python3

import numpy
import unittest
from afqmctools.analysis.extraction import extract_observable
from afqmctools.analysis.average import average_one_rdm


class TestRDM(unittest.TestCase):

    def test_average(self):
        f = 'qmc.s000.stat.h5'
        # Old format
        # base = 'Observables/BackPropagated/FullOneRDM/Average_1'
        rdm_av, rdm_errs = average_one_rdm(f, eqlb=1, ix=1)
        self.assertAlmostEqual(2*numpy.sum(rdm_av.real), 9.990239713872079)
        self.assertAlmostEqual(2*numpy.sum(rdm_av.imag), 0.009522316560114931)

    def test_extract(self):
        f = 'qmc.s000.stat.h5'
        # Old format
        # base = 'Observables/BackPropagated/FullOneRDM/Average_0'
        dm = 2*extract_observable(f, ix=0)
        self.assertAlmostEqual(numpy.sum(dm.real), 499.49878268185375)
        self.assertAlmostEqual(numpy.sum(dm.imag), 1.2453436020111393)

    def test_extract_single(self):
        f = 'qmc.s000.stat.h5'
        # base = 'Observables/BackPropagated/FullOneRDM/Average_3'
        dm = extract_observable(f, ix=3, sample=37)
        self.assertAlmostEqual(numpy.sum(dm).real, 5.00103256421043)


if __name__ == '__main__':
    import sys
    import os
    dir_path = os.path.dirname(os.path.realpath(__file__))
    sys.path.append(os.path.join(dir_path, '../'))
    from check_h1e_conv import plot_convergence
    plot_convergence('qmc.s000.stat.h5')
    unittest.main()
