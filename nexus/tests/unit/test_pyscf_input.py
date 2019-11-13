
import testing
from testing import value_eq,object_eq,text_eq


h2o_xyz = '''3

O  0.000000  0.000000  0.000000 
H  0.000000  0.757160  0.586260
H  0.000000  0.757160 -0.586260
'''


scf_template = '''#! /usr/bin/env python

from pyscf import scf

$system

mf = scf.RHF(mol)
mf.kernel()
'''


def test_import():
    from pyscf_input import PyscfInput,generate_pyscf_input
#end def test_import


def test_empty_init():
    from generic import obj
    from pyscf_input import PyscfInput,generate_pyscf_input

    ref = obj(
        addendum      = None,
        allow_not_set = set([]),
        checkpoint    = False,
        keywords      = set([]),
        prefix        = None,
        save_qmc      = False,
        template      = None,
        values        = obj(),
        )

    pi = PyscfInput()
    assert(object_eq(pi.to_obj(),ref))

    pi2 = generate_pyscf_input()
    assert(isinstance(pi2,PyscfInput))
    assert(object_eq(pi,pi2))
#end def test_empty_init



def test_generate():
    import os
    from generic import obj
    from physical_system import generate_physical_system
    from pyscf_input import generate_pyscf_input

    tpath = testing.setup_unit_test_output_directory('pyscf_input','test_generate')

    # water molecule
    xyz_path      = os.path.join(tpath,'H2O.xyz')
    template_path = os.path.join(tpath,'scf_template.py')

    open(xyz_path,'w').write(h2o_xyz)
    open(template_path,'w').write(scf_template)

    system = generate_physical_system(
        structure = xyz_path,
        )

    pi = generate_pyscf_input(
        template   = template_path,
        system     = system,
        mole       = obj(
            basis    = 'ccpvtz',
            symmetry = True,
            ),
        )

    ref_system = '''
        ### generated system text ###
        from pyscf import gto as gto_loc
        mol = gto_loc.Mole()
        mol.atom     = {0}
                       O    0.00000000   0.00000000   0.00000000
                       H    0.00000000   0.75716000   0.58626000
                       H    0.00000000   0.75716000  -0.58626000
                       {0}
        mol.basis    = 'ccpvtz'
        mol.unit     = 'A'
        mol.charge   = 0
        mol.spin     = 0
        mol.symmetry = True
        mol.build()
        ### end generated system text ###
        '''.format("'''")

    assert(pi.template is not None)
    assert(len(pi.values)==1 and 'system' in pi.values)
    assert(text_eq(pi.values.system,ref_system))

    ref_internal = obj(
        addendum      = None,
        allow_not_set = set([]),
        checkpoint    = False,
        keywords      = set(['system']),
        prefix        = None,
        save_qmc      = False,
        )
    del pi.template
    del pi.values
    assert(object_eq(pi.to_obj(),ref_internal))

    
    # diamond crystal
    system = generate_physical_system(
        units    = 'A',
        axes     = '''1.785   1.785   0.000
                      0.000   1.785   1.785
                      1.785   0.000   1.785''',
        elem_pos = '''
                   C  0.0000  0.0000  0.0000
                   C  0.8925  0.8925  0.8925
                   ''',
        kgrid    = (1,1,1),
        kshift   = (0,0,0),
        C        = 4,
        )


    pi = generate_pyscf_input(
        template   = template_path,
        system     = system,
        cell       = obj(
            basis         = 'bfd-vdz',
            ecp           = 'bfd',
            drop_exponent = 0.1,
            verbose       = 5,
            ),
        )

    ref_system = '''
        ### generated system text ###
        from numpy import array
        from pyscf.pbc import gto as gto_loc
        cell = gto_loc.Cell()
        cell.a             = {0}
                             1.78500000   1.78500000   0.00000000
                             0.00000000   1.78500000   1.78500000
                             1.78500000   0.00000000   1.78500000
                             {0}
        cell.basis         = 'bfd-vdz'
        cell.dimension     = 3
        cell.ecp           = 'bfd'
        cell.unit          = 'A'
        cell.atom          = {0}
                             C    0.00000000   0.00000000   0.00000000
                             C    0.89250000   0.89250000   0.89250000
                             {0}
        cell.drop_exponent = 0.1
        cell.verbose       = 5
        cell.charge        = 0
        cell.spin          = 0
        cell.build()
        kpts = array([
            [0.0, 0.0, 0.0]])
        ### end generated system text ###
        '''.format("'''")


    assert(pi.template is not None)
    assert(len(pi.values)==1 and 'system' in pi.values)
    assert(text_eq(pi.values.system,ref_system))

    del pi.template
    del pi.values
    assert(object_eq(pi.to_obj(),ref_internal))

#end def test_generate



def test_write():
    import os
    from generic import obj
    from physical_system import generate_physical_system
    from pyscf_input import generate_pyscf_input

    tpath = testing.setup_unit_test_output_directory('pyscf_input','test_write')

    # water molecule
    xyz_path      = os.path.join(tpath,'H2O.xyz')
    template_path = os.path.join(tpath,'scf_template.py')

    open(xyz_path,'w').write(h2o_xyz)
    open(template_path,'w').write(scf_template)

    system = generate_physical_system(
        structure = xyz_path,
        )

    pi = generate_pyscf_input(
        prefix     = 'scf',
        template   = template_path,
        system     = system,
        mole       = obj(
            verbose  = 5,
            basis    = 'ccpvtz',
            symmetry = True,
            ),
        save_qmc   = True,
        )

    write_path = os.path.join(tpath,'h2o.py')

    pi.write(write_path)

    assert(os.path.exists(write_path))

    text = open(write_path,'r').read()

    ref_text = '''
        #! /usr/bin/env python
        
        from pyscf import scf
        
        
        ### generated system text ###
        from pyscf import gto as gto_loc
        mol = gto_loc.Mole()
        mol.verbose  = 5
        mol.atom     = {0}
                       O    0.00000000   0.00000000   0.00000000
                       H    0.00000000   0.75716000   0.58626000
                       H    0.00000000   0.75716000  -0.58626000
                       {0}
        mol.basis    = 'ccpvtz'
        mol.unit     = 'A'
        mol.charge   = 0
        mol.spin     = 0
        mol.symmetry = True
        mol.build()
        ### end generated system text ###
        
        
        
        mf = scf.RHF(mol)
        mf.kernel()
        
        ### generated conversion text ###
        from PyscfToQmcpack import savetoqmcpack
        savetoqmcpack(mol,mf,'scf')
        ### end generated conversion text ###
        '''.format("'''")

    assert(text_eq(text,ref_text))


    # diamond crystal
    system = generate_physical_system(
        units    = 'A',
        axes     = '''1.785   1.785   0.000
                      0.000   1.785   1.785
                      1.785   0.000   1.785''',
        elem_pos = '''
                   C  0.0000  0.0000  0.0000
                   C  0.8925  0.8925  0.8925
                   ''',
        tiling   = (2,1,1),
        kgrid    = (1,1,1),
        kshift   = (0,0,0),
        C        = 4,
        )

    pi = generate_pyscf_input(
        prefix     = 'scf',
        template   = template_path,
        system     = system,
        cell       = obj(
            basis         = 'bfd-vdz',
            ecp           = 'bfd',
            drop_exponent = 0.1,
            verbose       = 5,
            ),
        save_qmc   = True,
        )

    write_path = os.path.join(tpath,'diamond.py')

    pi.write(write_path)

    assert(os.path.exists(write_path))

    text = open(write_path,'r').read()

    ref_text = '''
        #! /usr/bin/env python
        
        from pyscf import scf
        
        
        ### generated system text ###
        from numpy import array
        from pyscf.pbc import gto as gto_loc
        cell = gto_loc.Cell()
        cell.a             = {0}
                             1.78500000   1.78500000   0.00000000
                             0.00000000   1.78500000   1.78500000
                             1.78500000   0.00000000   1.78500000
                             {0}
        cell.basis         = 'bfd-vdz'
        cell.dimension     = 3
        cell.ecp           = 'bfd'
        cell.unit          = 'A'
        cell.atom          = {0}
                             C    0.00000000   0.00000000   0.00000000
                             C    0.89250000   0.89250000   0.89250000
                             {0}
        cell.drop_exponent = 0.1
        cell.verbose       = 5
        cell.charge        = 0
        cell.spin          = 0
        cell.build()
        kpts = array([
            [0.0, 0.0, 0.0] ,
            [0.4656748546088228, 0.4656748546088228, -0.4656748546088228]])
        ### end generated system text ###
        
        
        
        mf = scf.RHF(mol)
        mf.kernel()
        
        ### generated conversion text ###
        from PyscfToQmcpack import savetoqmcpack
        savetoqmcpack(cell,mf,'scf',kpts)
        ### end generated conversion text ###
        '''.format("'''")

    text = text.replace('[',' [ ').replace(']',' ] ')
    ref_text = ref_text.replace('[',' [ ').replace(']',' ] ')

    assert(text_eq(text,ref_text))

#end def test_write
