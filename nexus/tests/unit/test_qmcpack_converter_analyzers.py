


def test_import():
    from qmcpack_converters import Pw2qmcpackAnalyzer
    from qmcpack_converters import Convert4qmcAnalyzer
    from qmcpack_converters import PyscfToAfqmcAnalyzer
#end def import



def test_pw2qmcpack_analyzer_init():
    from qmcpack_converters import Pw2qmcpackAnalyzer
    # empty init
    pa = Pw2qmcpackAnalyzer(None)
#end def test_pw2qmcpack_analyzer_init



def test_convert4qmc_analyzer_init():
    from qmcpack_converters import Convert4qmcAnalyzer
    # empty init
    pa = Convert4qmcAnalyzer(None)
#end def test_convert4qmc_analyzer_init



def test_pyscf_to_afqmc_analyzer_init():
    from qmcpack_converters import PyscfToAfqmcAnalyzer
    # empty init
    pa = PyscfToAfqmcAnalyzer(None)
#end def test_pyscf_to_afqmc_analyzer_init
