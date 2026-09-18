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

    analyzer = analyze_output(code,outfile=filepath)
    assert(analyzer.__class__.__name__==analyzer_name)
#end def test_analyze_output_constructs_analyzer


def test_analyze_output_unified_paths():
    from .. import analyze_output

    pwscf_path = TEST_DIR/'test_pwscf_analyzer_files/qe_7_0/default/scf'
    pwscf = analyze_output(
        'pwscf',
        input   = 'pwscf.in',
        outfile = 'pwscf.out',
        path    = pwscf_path,
        )
    assert pwscf.analysis_state=='analyzed'
    assert pwscf.input is not None

    rmg_path = TEST_DIR/'test_rmg_analyzer_files/electronic'
    rmg = analyze_output(
        'rmg',
        input   = 'input.scf',
        outfile = 'input.scf.02.log',
        path    = rmg_path,
        )
    assert rmg.analysis_state=='analyzed'
    assert rmg.input is not None

    deferred = analyze_output('pwscf',analyze=False,strict=False)
    assert deferred.analysis_state=='not_analyzed'
#end def test_analyze_output_unified_paths
