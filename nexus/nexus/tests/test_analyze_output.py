import pytest

from . import TEST_DIR, NexusTestOrder


pytestmark = pytest.mark.order(NexusTestOrder.ANALYZE_OUTPUT)


@pytest.mark.parametrize(
    'code,filepath,analyzer_name',
    (
        (
            'pwscf',
            TEST_DIR/'test_pwscf_analyzer_files/qe_7_0/default/scf/pwscf.out',
            'PwscfAnalyzer',
            ),
        (
            'rmg',
            TEST_DIR/'test_rmg_analyzer_files/electronic/input.scf.02.log',
            'RmgAnalyzer',
            ),
        ),
    )
def test_analyze_output_constructs_analyzer(code,filepath,analyzer_name):
    from .. import analyze_output

    analyzer = analyze_output(code,path=filepath)
    assert(analyzer.__class__.__name__==analyzer_name)
#end def test_analyze_output_constructs_analyzer
