
import testing
from testing import value_eq,object_eq


def test_import():
    from qmcpack_converters import Pw2qmcpackAnalyzer
    from qmcpack_converters import Pw2qmcpack
    from qmcpack_converters import generate_pw2qmcpack

    from qmcpack_converters import Convert4qmcAnalyzer
    from qmcpack_converters import Convert4qmc
    from qmcpack_converters import generate_convert4qmc

    from qmcpack_converters import PyscfToAfqmcAnalyzer

#end def test_import




def test_pyscf_to_afqmc_analyzer_init():
    from qmcpack_converters import PyscfToAfqmcAnalyzer
    # empty init
    pa = PyscfToAfqmcAnalyzer(None)
#end def test_pyscf_to_afqmc_analyzer_init

