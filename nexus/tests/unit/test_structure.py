#!/env/bin/python

import numpy as np
from testing import value_eq

def test_count_kshells():
    from test_physical_system import example_structure_h4
    s1 = example_structure_h4()
    kf = 1.465
    kcut = 5*kf
    nksh = s1.count_kshells(kcut)
    assert nksh == 13
#end def test_count_kshells


def test_rinscribe():
    from structure import generate_structure
    s = generate_structure(
        structure = 'graphene',
        cell      = 'prim',
        tiling    = (4,4,1),
        )
    ri = s.rinscribe()
    assert(value_eq(float(ri),4.2643090882345795))
#end def test_rinscribe


if __name__ == '__main__':
  test_count_kshells()
# end __main__
