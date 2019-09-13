import numpy
import unittest

from afqmctools.analysis.rdm import (
        get_one_rdm_av,
        extract_rdm,
        extract_rdm_single
        )


class TestRDM(unittest.TestCase):

    def test_average(self):
        f = 'qmc.s000.scalar.h5'
        base = 'Observables/BackPropagated/FullOneRDM/Average_1'
        rdm_av, rdm_errs = get_one_rdm_av(f, 1, dm_name=base)
        self.assertAlmostEqual(numpy.sum(rdm_av.real), 9.990239713872079)
        self.assertAlmostEqual(numpy.sum(rdm_av.imag), 0.009522316560114931)

    def test_extract(self):
        f = 'qmc.s000.scalar.h5'
        base = 'Observables/BackPropagated/FullOneRDM/Average_0'
        dm, weights = extract_rdm(f, dm_name=base)
        self.assertAlmostEqual(numpy.sum(dm.real), 499.49878268185375)
        self.assertAlmostEqual(numpy.sum(dm.imag), 1.2453436020111393)

    def test_extract_single(self):
        f = 'qmc.s000.scalar.h5'
        base = 'Observables/BackPropagated/FullOneRDM/Average_3'
        dm, weight = extract_rdm_single(f, 37, dm_name=base)
        self.assertAlmostEqual(numpy.sum(dm).real, 160)


if __name__ == '__main__':
    from check_h1e_conv import plot_convergence
    plot_convergence('qmc.s000.scalar.h5')
    unittest.main()
