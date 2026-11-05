import pytest

from ..qmcpack_converters import Convert4qmcAnalyzer, Pw2qmcpackAnalyzer, PyscfToAfqmcAnalyzer
from . import NexusTestOrder

pytestmark = pytest.mark.order(NexusTestOrder.QMCPACK_CONVERTER_ANALYZERS)

def test_pw2qmcpack_analyzer_init():
    # empty init
    pa = Pw2qmcpackAnalyzer(None)
#end def test_pw2qmcpack_analyzer_init



def test_convert4qmc_analyzer_init():
    # empty init
    pa = Convert4qmcAnalyzer(None)
#end def test_convert4qmc_analyzer_init



def test_pyscf_to_afqmc_analyzer_init():
    # empty init
    pa = PyscfToAfqmcAnalyzer(None)
#end def test_pyscf_to_afqmc_analyzer_init
