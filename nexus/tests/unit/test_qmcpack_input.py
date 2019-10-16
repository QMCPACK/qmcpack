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


def test_excited_state():
    from nexus import generate_physical_system
    from nexus import generate_qmcpack_input

    dia = generate_physical_system(
        units     = 'A',
        axes      = [[ 1.785,  1.785,  0.   ],
                     [ 0.   ,  1.785,  1.785],
                     [ 1.785,  0.   ,  1.785]],
        elem      = ['C','C'],
        pos       = [[ 0.    ,  0.    ,  0.    ],
                     [ 0.8925,  0.8925,  0.8925]],
        tiling    = [3,1,3], 
        kgrid     = (1,1,1), 
        kshift    = (0,0,0), 
        C         = 4
        )
  
  
    # test kp_index, band_index format (format="band")
    qmc_optical = generate_qmcpack_input(
        det_format     = 'old',
        spin_polarized = True,
        system         = dia,
        excitation     = ['up', '0 3 4 4'], #
        pseudos        = ['C.BFD.xml'],
        jastrows       = [],
        qmc            = 'vmc',
        )

    expect = '''<slaterdeterminant>
   <determinant id="updet" size="36">
      <occupation mode="excited" spindataset="0" pairs="1" format="band">             
0 3 4 4
       </occupation>
   </determinant>
   <determinant id="downdet" size="36">
      <occupation mode="ground" spindataset="1"/>
   </determinant>
</slaterdeterminant>'''.strip()

    text = qmc_optical.get('slaterdeterminant').write().strip()
    assert(text==expect)


    # test energy_index (format="energy")
    qmc_optical = generate_qmcpack_input(
        det_format     = 'old',
        spin_polarized = True,
        system         = dia,
        excitation     = ['up', '-35 36'], #
        pseudos        = ['C.BFD.xml'],
        jastrows       = [],
        qmc            = 'vmc',
        )

    expect = '''<slaterdeterminant>
   <determinant id="updet" size="36">
      <occupation mode="excited" spindataset="0" pairs="1" format="energy">             
-35 36
        </occupation>
   </determinant>
   <determinant id="downdet" size="36">
      <occupation mode="ground" spindataset="1"/>
   </determinant>
</slaterdeterminant>'''.strip()

    text = qmc_optical.get('slaterdeterminant').write().strip()
    assert(text==expect)

#end def test_excited_state


def test_symbolic_excited_state():
    from nexus import generate_physical_system
    from nexus import generate_qmcpack_input

    # remove once explicit checks/guards for seekpath are made
    return

    dia = generate_physical_system(
        units     = 'A',
        axes      = [[ 1.785,  1.785,  0.   ],
                     [ 0.   ,  1.785,  1.785],
                     [ 1.785,  0.   ,  1.785]],
        elem      = ['C','C'],
        pos       = [[ 0.    ,  0.    ,  0.    ],
                     [ 0.8925,  0.8925,  0.8925]],
        use_prim  = True,    # Use SeeK-path library to identify prim cell
        tiling    = [2,1,2], 
        kgrid     = (1,1,1), 
        kshift    = (0,0,0), # Assumes we study transitions from Gamma. For non-gamma tilings, use kshift appropriately
        #C         = 4
        )

  
    qmc_optical = generate_qmcpack_input(
        det_format     = 'old',
        input_type     = 'basic',
        spin_polarized = True,
        system         = dia,
        excitation     = ['up', 'gamma vb x cb'], 
        jastrows       = [],
        qmc            = 'vmc',
        )

    expect = '''<slaterdeterminant>
   <determinant id="updet" size="24">
      <occupation mode="excited" spindataset="0" pairs="1" format="band">             
0 5 3 6
       </occupation>
   </determinant>
   <determinant id="downdet" size="24">
      <occupation mode="ground" spindataset="1"/>
   </determinant>
</slaterdeterminant>'''.strip()
    text = qmc_optical.get('slaterdeterminant').write().strip()
    assert(text==expect)
  

    qmc_optical = generate_qmcpack_input(
        det_format     = 'old',
        input_type     = 'basic',
        spin_polarized = True,
        system         = dia,
        excitation     = ['up', 'gamma vb-1 x cb'], 
        jastrows       = [],
        qmc            = 'vmc',
        )

    expect = '''<slaterdeterminant>
   <determinant id="updet" size="24">
      <occupation mode="excited" spindataset="0" pairs="1" format="band">             
0 4 3 6
       </occupation>
   </determinant>
   <determinant id="downdet" size="24">
      <occupation mode="ground" spindataset="1"/>
   </determinant>
</slaterdeterminant>'''.strip()
    text = qmc_optical.get('slaterdeterminant').write().strip()
    assert(text==expect)


    qmc_optical = generate_qmcpack_input(
        det_format     = 'old',
        input_type     = 'basic',
        spin_polarized = True,
        system         = dia,
        excitation     = ['up', 'gamma vb x cb+1'], 
        jastrows       = [],
        qmc            = 'vmc',
        )

    expect = '''<slaterdeterminant>
   <determinant id="updet" size="24">
      <occupation mode="excited" spindataset="0" pairs="1" format="band">             
0 5 3 7
       </occupation>
   </determinant>
   <determinant id="downdet" size="24">
      <occupation mode="ground" spindataset="1"/>
   </determinant>
</slaterdeterminant>'''.strip()
    text = qmc_optical.get('slaterdeterminant').write().strip()
    assert(text==expect)

#end def test_symbolic_excited_state



if __name__ == '__main__':
    test_generate_kspace_jastrow()
#end __main__
