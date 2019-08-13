#!/usr/bin/env python
def test_generate_kspace_jastrow():
  from qmcpack_input import generate_kspace_jastrow
  kjas = generate_kspace_jastrow(1.0, 2.0, 2, 4)
  expect = '''<jastrow type="kSpace" name="Jk" source="ion0">
   <correlation kc="1.0" type="One-Body" symmetry="isotropic">
      <coefficients id="cG1" type="Array">         
0 0
      </coefficients>
   </correlation>
   <correlation kc="2.0" type="Two-Body" symmetry="isotropic">
      <coefficients id="cG2" type="Array">         
0 0 0 0
      </coefficients>
   </correlation>
</jastrow>
'''
  text = kjas.write()
  assert text == expect
#end def test_generate_kspace_jastrow

if __name__ == '__main__':
  test_generate_kspace_jastrow()
#end __main__
