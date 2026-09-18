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


def test_analyze_output_simulation_forms(tmp_path):
    from .. import analyze_output
    from ..machines import job
    from ..pwscf import generate_pwscf
    from ..testing import clear_all_sims
    from .test_pwscf_simulation import get_system

    sim = generate_pwscf(
        job    = job(machine='ws1',cores=1),
        system = get_system(),
        )
    sentinel = object()
    sim.imresdir = str(tmp_path)
    imagepath = tmp_path/sim.analyzer_image
    imagepath.touch()
    sim.load_analyzer_image = lambda:sentinel

    assert analyze_output(sim) is sentinel
    assert analyze_output(input=sim) is sentinel
    assert analyze_output('pwscf',sim) is sentinel
    with pytest.raises(ValueError,match='does not match'):
        analyze_output('rmg',sim)
    with pytest.raises(ValueError,match='additional arguments'):
        analyze_output(sim,analyze=False)

    clear_all_sims()
#end def test_analyze_output_simulation_forms
