import pytest
from . import NexusTestOrder
pytestmark = pytest.mark.order(NexusTestOrder.QMC_FIT)

from ..generic import generic_settings
generic_settings.raise_error = True

from . import TEST_DIR
from .. import versions
from ..testing import execute,text_eq



if versions.scipy_available:
    def test_fit(tmp_path):

        exe = TEST_DIR.parent / "bin/qmc-fit"        
        dmc_path = TEST_DIR / "test_qmcpack_analyzer_files/diamond_gamma/dmc"

        dmc_infile = dmc_path / 'dmc.in.xml'
        assert(dmc_infile.exists())

        command = f"{exe} ts --noplot -e 10 -s 1 -t '0.02 0.01 0.005' -f linear {dmc_path}/*scalar*"

        out,err,rc = execute(command)

        out_ref = '''
            fit function  : linear
            fitted formula: (-10.5271 +/- 0.0021) + (-0.28 +/- 0.17)*t
            intercept     : -10.5271 +/- 0.0021  Ha
            '''
        
        def process_text(t):
            return t.replace('(',' ( ').replace(')',' ) ')
        #end def process_text

        out = process_text(out)
        out_ref = process_text(out_ref)

        assert(text_eq(out,out_ref,atol=1e-2,rtol=1e-2))

    #end def test_fit
#end if
