
import testing
from testing import value_eq,object_eq,text_eq



def test_import():
    from qmcpack_converters import Pw2qmcpackInput
    from qmcpack_converters import generate_pw2qmcpack_input

    from qmcpack_converters import Convert4qmcInput
    from qmcpack_converters import generate_convert4qmc_input

    from qmcpack_converters import PyscfToAfqmcInput
    from qmcpack_converters import generate_pyscf_to_afqmc_input
#end def test_import




pw2qmcpack_in = '''&inputpp
  prefix = 'pwscf'
  outdir = 'pwscf_output'
/
'''

def test_pw2qmcpack_empty_init():
    from qmcpack_converters import Pw2qmcpackInput
    from qmcpack_converters import generate_pw2qmcpack_input

    pi = Pw2qmcpackInput()

    pi = generate_pw2qmcpack_input()
    assert(isinstance(pi,Pw2qmcpackInput))
#end def test_pw2qmcpack_empty_init



def test_pw2qmcpack_read():
    import os
    from generic import obj
    from qmcpack_converters import Pw2qmcpackInput

    tpath = testing.setup_unit_test_output_directory('qmcpack_converter_input','test_pw2qmcpack_read')

    infile_path = os.path.join(tpath,'p2q.in')
    open(infile_path,'w').write(pw2qmcpack_in)

    pi = Pw2qmcpackInput(infile_path)
    
    pi_ref = obj(
        inputpp = obj(
            prefix     = 'pwscf',
            outdir     = 'pwscf_output',
            ),
        )

    assert(object_eq(pi.to_obj(),pi_ref))
#end def test_pw2qmcpack_read



def test_pw2qmcpack_write():
    import os
    from generic import obj
    from qmcpack_converters import Pw2qmcpackInput

    tpath = testing.setup_unit_test_output_directory('qmcpack_converter_input','test_pw2qmcpack_write')

    infile_path = os.path.join(tpath,'p2q.in')
    open(infile_path,'w').write(pw2qmcpack_in)

    write_path = os.path.join(tpath,'p2q_write.in')
    pi_write = Pw2qmcpackInput(infile_path)
    
    pi_write.write(write_path)

    pi_read = Pw2qmcpackInput(write_path)

    pi_ref = obj(
        inputpp = obj(
            prefix = 'pwscf',
            outdir = 'pwscf_output',
            ),
        )

    assert(object_eq(pi_read.to_obj(),pi_ref))
#end def test_pw2qmcpack_write



def test_pw2qmcpack_generate():
    from generic import obj
    from qmcpack_converters import generate_pw2qmcpack_input

    pi = generate_pw2qmcpack_input(
        prefix = 'pwscf',
        outdir = 'pwscf_output',
        )

    pi_ref = obj(
        inputpp = obj(
            prefix     = 'pwscf',
            outdir     = 'pwscf_output',
            write_psir = False,
            ),
        )

    assert(object_eq(pi.to_obj(),pi_ref))
#end def test_pw2qmcpack_generate



def test_pyscf_to_afqmc_input_init():
    from qmcpack_converters import PyscfToAfqmcInput
    from qmcpack_converters import generate_pyscf_to_afqmc_input

    # empty init
    pi = PyscfToAfqmcInput()
    assert(pi.is_valid())
    assert(set(pi.keys())==set(PyscfToAfqmcInput.input_defaults.keys()))
    for v in pi:
        assert(v is None or v==False or v=='pyscf_to_afqmc.py')
    #end for

    pi2 = generate_pyscf_to_afqmc_input()
    assert(pi2.is_valid())
    assert(object_eq(pi2,pi))

    # simple init
    pi = generate_pyscf_to_afqmc_input(
        input              = 'scf.chk',
        output             = 'afqmc.h5',
        cholesky_threshold = 1e-5,
        verbose            = True,
        )
    assert(pi.is_valid())
    assert(pi.input=='scf.chk')
    assert(pi.output=='afqmc.h5')
    assert(value_eq(pi.cholesky_threshold,1e-5))
    assert(value_eq(pi.verbose,True))

    pi2 = generate_pyscf_to_afqmc_input(
        i = 'scf.chk',
        o = 'afqmc.h5',
        t = 1e-5,
        v = True,
        )
    assert(pi2.is_valid())
    assert(object_eq(pi2,pi))

#end def test_pyscf_to_afqmc_input_init



def test_pyscf_to_afqmc_input_write():
    from qmcpack_converters import generate_pyscf_to_afqmc_input

    # example 1: neon atom
    pi = generate_pyscf_to_afqmc_input(
        i = 'scf.chk',
        o = 'afqmc.h5',
        t = 1e-5,
        v = True,
        )
    ref = 'pyscf_to_afqmc.py -i scf.chk -o afqmc.h5 -t 1e-05 -v'
    assert(pi.write()==ref)

    # example 2: frozen core
    pi = generate_pyscf_to_afqmc_input(
        i = 'scf.chk',
        o = 'afqmc.h5',
        t = 1e-5,
        c = (8,22),
        v = True,
        )
    ref = 'pyscf_to_afqmc.py -i scf.chk -o afqmc.h5 -t 1e-05 -c 8,22 -v'
    assert(pi.write()==ref)

    # example 3: uhf trial
    pi = generate_pyscf_to_afqmc_input(
        i = 'scf.chk',
        o = 'afqmc.h5',
        t = 1e-5,
        a = True,
        v = True,
        )
    ref = 'pyscf_to_afqmc.py -i scf.chk -o afqmc.h5 -t 1e-05 -a -v'
    assert(pi.write()==ref)

    # example 8: diamond k-point symmetry
    pi = generate_pyscf_to_afqmc_input(
        i = 'scf.chk',
        o = 'afqmc.h5',
        t = 1e-5,
        a = True,
        k = True,
        v = True,
        )
    ref = 'pyscf_to_afqmc.py -i scf.chk -o afqmc.h5 -t 1e-05 -k -a -v'
    assert(pi.write()==ref)

#end def test_pyscf_to_afqmc_input_write
