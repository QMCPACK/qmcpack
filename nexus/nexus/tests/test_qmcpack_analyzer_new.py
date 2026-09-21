import numpy as np
import pytest

from . import NexusTestOrder, TEST_DIR


pytestmark = pytest.mark.order(NexusTestOrder.QMCPACK_ANALYZER)


ANALYZER_FILES = TEST_DIR / 'test_qmcpack_analyzer_files'
INPUT_FILES = TEST_DIR / 'test_qmcpack_input_files'


def test_scalar_info():
    from ..qmcpack_analyzer_new import scalar_info

    assert scalar_info.aliases.LocalEnergy == 'E'
    assert scalar_info.inv_aliases.E == 'LocalEnergy'
    assert scalar_info.integer == {'TotalSamples', 'NumOfWalkers'}
    assert scalar_info.constant == {'IonIon', 'KEcorr', 'MPC'}
    assert scalar_info.nonenergy <= scalar_info.analyze


@pytest.mark.parametrize(
    'relative_path,nrows',
    [
        ('diamond_gamma/vmc/vmc.s000.scalar.dat', 100),
        ('diamond_gamma/dmc/dmc.s000.scalar.dat', 10),
        ('diamond_gamma/opt/opt.s000.scalar.dat', 50),
        ('diamond_twist/vmc/vmc.g002.s000.scalar.dat', 100),
    ],
)
def test_read_scalar_file(relative_path, nrows):
    from ..qmcpack_analyzer_new import read_scalar_file

    filepath = ANALYZER_FILES / relative_path
    data, issues = read_scalar_file(
        str(filepath), issues=True, add_variance=True, remove_index=True
    )

    assert not issues.failed()
    assert issues.issue_set() == set()
    assert 'index' not in data
    assert 'LocalEnergy_sq' not in data
    assert 'Variance' in data
    assert len(data) == 12
    assert {len(values) for values in data.values()} == {nrows}
    assert all(np.isfinite(values).all() for values in data.values())
    raw_data = read_scalar_file(str(filepath))
    np.testing.assert_allclose(
        data['Variance'], raw_data['LocalEnergy_sq'] - raw_data['LocalEnergy']**2
    )


def test_read_scalar_file_missing():
    from ..qmcpack_analyzer_new import read_scalar_file

    data, issues = read_scalar_file(
        str(ANALYZER_FILES / 'missing.scalar.dat'), issues=True
    )

    assert data == {}
    assert issues.issue_set() == {'no_file'}
    assert issues.failed()


@pytest.mark.parametrize(
    'qmc,prefix,series,group_index,expected',
    [
        ('vmc', 'vmc', 0, None, ('vmc.s000.scalar.dat',)),
        ('dmc', 'dmc', 3, None,
         ('dmc.s003.scalar.dat', 'dmc.s003.dmc.dat')),
        ('opt', 'opt', 4, None,
         ('opt.s004.scalar.dat', 'opt.s004.opt.xml', 'opt.s004.vp.h5')),
        ('vmc', 'vmc', 0, 2, ('vmc.g002.s000.scalar.dat',)),
    ],
)
def test_qmcpack_analyzer_outfiles(qmc, prefix, series, group_index, expected):
    from ..qmcpack_analyzer_new import qmcpack_analyzer_outfiles

    assert qmcpack_analyzer_outfiles(qmc, prefix, series, group_index) == expected


@pytest.mark.parametrize(
    'relative_path,qmc_type,prefix,group_index,nseries,expected_outfiles',
    [
        ('diamond_gamma/vmc/vmc.in.xml', 'vmc', 'vmc', None, 1,
         ('vmc.s000.scalar.dat',)),
        ('diamond_gamma/dmc/dmc.in.xml', 'dmc', 'dmc', None, 4,
         ('dmc.s003.scalar.dat', 'dmc.s003.dmc.dat')),
        ('diamond_gamma/opt/opt.in.xml', 'opt', 'opt', None, 6,
         ('opt.s005.scalar.dat', 'opt.s005.opt.xml', 'opt.s005.vp.h5')),
        ('diamond_twist/vmc/vmc.g002.twistnum_2.in.xml', 'vmc', 'vmc', 2, 1,
         ('vmc.g002.s000.scalar.dat',)),
    ],
)
def test_qmcpack_input_info_from_run_inputs(
    relative_path, qmc_type, prefix, group_index, nseries, expected_outfiles
):
    from ..qmcpack_analyzer_new import QmcpackInputInfo

    info = QmcpackInputInfo(str(ANALYZER_FILES / relative_path))

    assert info.qmc_type == qmc_type
    assert info.prefix == prefix
    assert info.group_index == group_index
    assert info.has_twist
    assert len(info.qmc_info) == nseries
    assert info.qmc_info[nseries - 1].outfiles == expected_outfiles


def test_qmcpack_input_info_from_input_fixture():
    from ..qmcpack_analyzer_new import QmcpackInputInfo

    info = QmcpackInputInfo(str(INPUT_FILES / 'VO2_M1_afm.in.xml'))

    assert info.qmc_type == 'dmc'
    assert info.prefix == 'qmc'
    assert info.group_index is None
    assert info.has_twist
    assert list(info.qmc_info) == [0, 1, 2]
    assert info.qmc_info[0].qmc == 'vmc'
    assert info.qmc_info[2].outfiles == (
        'qmc.s002.scalar.dat', 'qmc.s002.dmc.dat'
    )
