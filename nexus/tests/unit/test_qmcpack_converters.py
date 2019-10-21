
import testing
from testing import value_eq,object_eq


def test_import():
    from qmcpack_converters import Pw2qmcpackInput
    from qmcpack_converters import Pw2qmcpackAnalyzer
    from qmcpack_converters import Pw2qmcpack
    from qmcpack_converters import generate_pw2qmcpack_input
    from qmcpack_converters import generate_pw2qmcpack

    from qmcpack_converters import Convert4qmcInput
    from qmcpack_converters import Convert4qmcAnalyzer
    from qmcpack_converters import Convert4qmc
    from qmcpack_converters import generate_convert4qmc_input
    from qmcpack_converters import generate_convert4qmc

    from qmcpack_converters import PyscfToAfqmcInput
    from qmcpack_converters import PyscfToAfqmcAnalyzer
    from qmcpack_converters import generate_pyscf_to_afqmc_input

#end def test_import



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



def test_pyscf_to_afqmc_analyzer_init():
    from qmcpack_converters import PyscfToAfqmcAnalyzer
    # empty init
    pa = PyscfToAfqmcAnalyzer(None)
#end def test_pyscf_to_afqmc_analyzer_init

