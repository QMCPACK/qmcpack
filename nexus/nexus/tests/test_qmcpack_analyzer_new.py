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


def test_read_scalar_file_single_column(tmp_path):
    from ..qmcpack_analyzer_new import read_scalar_file

    one_row = tmp_path / 'one_row.scalar.dat'
    one_row.write_text('# LocalEnergy\n-1.25\n')
    data, issues = read_scalar_file(str(one_row), issues=True)

    assert issues.issue_set() == set()
    np.testing.assert_array_equal(data['LocalEnergy'], [-1.25])

    multiple_rows = tmp_path / 'multiple_rows.scalar.dat'
    multiple_rows.write_text('# LocalEnergy\n-1.0\n-2.0\n-3.0\n')
    data, issues = read_scalar_file(str(multiple_rows), issues=True)

    assert issues.issue_set() == set()
    np.testing.assert_array_equal(data['LocalEnergy'], [-1.0, -2.0, -3.0])


@pytest.mark.parametrize(
    'header,rows',
    [
        ('# index LocalEnergy Kinetic\n', '0 -1.0\n1 -2.0\n'),
        ('# index LocalEnergy\n', '0 -1.0 2.0\n1 -2.0 3.0\n'),
    ],
)
def test_read_scalar_file_bad_column_count(tmp_path, header, rows):
    from ..qmcpack_analyzer_new import read_scalar_file

    filepath = tmp_path / 'bad_col_count.scalar.dat'
    filepath.write_text(header + rows)
    data, issues = read_scalar_file(str(filepath), issues=True)

    assert data == {}
    assert issues.issue_set() == {'bad_col_count', 'no_usable_vals'}
    assert issues.failed()


def test_read_scalar_file_nan_rows_and_trailing_columns(tmp_path):
    from ..qmcpack_analyzer_new import read_scalar_file

    filepath = tmp_path / 'nan_rows.scalar.dat'
    filepath.write_text(
        '# index LocalEnergy Kinetic ignored LocalPotential\n'
        '0   -1.0  2.0  nan  3.0\n'
        '1    nan  4.0  5.0  6.0\n'
        '2   -3.0  6.0  7.0  8.0\n'
    )

    data, issues = read_scalar_file(str(filepath), issues=True)
    assert set(data) == {'index', 'LocalEnergy', 'Kinetic'}
    assert issues.issue_set() == {'nan_vals'}
    np.testing.assert_array_equal(data['index'], [0.0, 2.0])
    np.testing.assert_array_equal(data['LocalEnergy'], [-1.0, -3.0])

    data, issues = read_scalar_file(
        str(filepath), issues=True, trim_nan=False
    )
    assert issues.issue_set() == {'nan_vals'}
    assert len(data['LocalEnergy']) == 3
    assert np.isnan(data['LocalEnergy'][1])


def test_read_scalar_file_ignores_nan_outside_scalar_columns(tmp_path):
    from ..qmcpack_analyzer_new import read_scalar_file

    filepath = tmp_path / 'ignored_nan.scalar.dat'
    filepath.write_text(
        '# index LocalEnergy ignored\n'
        'nan -1.0 nan\n'
    )

    data, issues = read_scalar_file(str(filepath), issues=True)

    assert issues.issue_set() == set()
    assert set(data) == {'index', 'LocalEnergy'}
    assert np.isnan(data['index'][0])
    np.testing.assert_array_equal(data['LocalEnergy'], [-1.0])


def test_read_scalar_file_corrupt_end(tmp_path):
    from ..qmcpack_analyzer_new import read_scalar_file

    filepath = tmp_path / 'corrupt_end.scalar.dat'
    filepath.write_text(
        '# index LocalEnergy Kinetic\n'
        '0 -1.0 2.0\n'
        '1 -2.0 3.0\n'
        '2 unfinished nonsense\n'
        'additional trailing garbage\n'
    )

    data, issues = read_scalar_file(str(filepath), issues=True, nrows=2)

    assert issues.issue_set() == {'corrupt_end'}
    assert not issues.failed()
    np.testing.assert_array_equal(data['index'], [0.0, 1.0])


def test_read_scalar_file_middle_corruption(tmp_path):
    from ..qmcpack_analyzer_new import read_scalar_file

    filepath = tmp_path / 'middle_corruption.scalar.dat'
    filepath.write_text(
        '# index LocalEnergy Kinetic\n'
        '0 -1.0 2.0\n'
        '1 nonsense 3.0\n'
        '2 -3.0 4.0\n'
    )

    data, issues = read_scalar_file(str(filepath), issues=True, nrows=2)

    assert issues.issue_set() == {'unparsable_vals'}
    assert issues.failed()
    np.testing.assert_array_equal(data['index'], [0.0])


def test_read_scalar_file_nrows_uses_complete_numeric_rows(tmp_path):
    from ..qmcpack_analyzer_new import read_scalar_file

    filepath = tmp_path / 'incomplete.scalar.dat'
    filepath.write_text(
        '# index LocalEnergy Kinetic\n'
        '0 -1.0 2.0\n'
        '1 -2.0 3.0\n'
        'garbage that makes the physical line count long enough\n'
    )

    data, issues = read_scalar_file(str(filepath), issues=True, nrows=3)

    assert len(data['index']) == 2
    assert issues.issue_set() == {'corrupt_end', 'incomplete'}
    assert not issues.failed()

    nan_file = tmp_path / 'nan_count.scalar.dat'
    nan_file.write_text(
        '# index LocalEnergy Kinetic\n'
        '0 -1.0 2.0\n'
        '1 nan 3.0\n'
    )
    data, issues = read_scalar_file(str(nan_file), issues=True, nrows=2)

    assert len(data['index']) == 1
    assert issues.issue_set() == {'nan_vals'}

    data, issues = read_scalar_file(str(nan_file), issues=True, nrows=1)
    assert issues.issue_set() == {'nan_vals', 'incomplete'}


def test_read_scalar_file_final_short_numeric_row(tmp_path):
    from ..qmcpack_analyzer_new import read_scalar_file

    filepath = tmp_path / 'short_final.scalar.dat'
    filepath.write_text(
        '# index LocalEnergy Kinetic\n'
        '0 -1.0 2.0\n'
        '1 -2.0\n'
    )

    data, issues = read_scalar_file(str(filepath), issues=True, nrows=2)

    assert issues.issue_set() == {'incomplete'}
    assert not issues.failed()
    np.testing.assert_array_equal(data['index'], [0.0])

    filepath.write_text(
        '# index LocalEnergy Kinetic\n'
        '0 -1.0 2.0\n'
        '1 -2.0\n'
        'trailing garbage with a different width\n'
    )
    data, issues = read_scalar_file(str(filepath), issues=True, nrows=2)

    assert issues.issue_set() == {'corrupt_end', 'incomplete'}
    assert not issues.failed()
    np.testing.assert_array_equal(data['index'], [0.0])

    filepath.write_text(
        '# index LocalEnergy Kinetic\n'
        '0 -1.0 2.0\n'
        'garbage before the partial row\n'
        '1 -2.0\n'
    )
    data, issues = read_scalar_file(str(filepath), issues=True, nrows=2)

    assert issues.issue_set() == {'corrupt_end', 'incomplete'}
    assert not issues.failed()
    np.testing.assert_array_equal(data['index'], [0.0])


def test_read_scalar_file_bad_nonfinal_widths(tmp_path):
    from ..qmcpack_analyzer_new import read_scalar_file

    short_middle = tmp_path / 'short_middle.scalar.dat'
    short_middle.write_text(
        '# index LocalEnergy Kinetic\n'
        '0 -1.0 2.0\n'
        '1 -2.0\n'
        '2 -3.0 4.0\n'
    )
    data, issues = read_scalar_file(str(short_middle), issues=True, nrows=2)

    assert issues.issue_set() == {'uneven_cols'}
    assert issues.failed()
    np.testing.assert_array_equal(data['index'], [0.0])

    long_final = tmp_path / 'long_final.scalar.dat'
    long_final.write_text(
        '# index LocalEnergy Kinetic\n'
        '0 -1.0 2.0\n'
        '1 -2.0 3.0 4.0\n'
    )
    data, issues = read_scalar_file(str(long_final), issues=True, nrows=2)

    assert issues.issue_set() == {'uneven_cols', 'incomplete'}
    assert issues.failed()
    np.testing.assert_array_equal(data['index'], [0.0])


def test_read_scalar_file_unparsable_rows_affect_width_check(tmp_path):
    from ..qmcpack_analyzer_new import read_scalar_file

    filepath = tmp_path / 'unparsable_width.scalar.dat'
    filepath.write_text(
        '# index LocalEnergy Kinetic\n'
        '0 -1.0 2.0\n'
        'nonsense\n'
        '2 -3.0 4.0\n'
    )

    data, issues = read_scalar_file(str(filepath), issues=True)

    assert issues.issue_set() == {'unparsable_vals', 'uneven_cols'}
    assert issues.failed()
    np.testing.assert_array_equal(data['index'], [0.0])


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


def test_qmcpack_input_info_method_aliases_and_multi_qmc_loop(tmp_path):
    from ..qmcpack_analyzer_new import QmcpackInputInfo

    filepath = tmp_path / 'methods.in.xml'
    filepath.write_text(
        '<simulation>\n'
        '  <project id="methods" series="2" />\n'
        '  <loop max="2">\n'
        '    <qmc method="cslinear">\n'
        '      <parameter name="blocks">3</parameter>\n'
        '    </qmc>\n'
        '    <qmc method="vmc_batch">\n'
        '      <parameter name="blocks">4</parameter>\n'
        '    </qmc>\n'
        '  </loop>\n'
        '  <qmc method="dmc_batch">\n'
        '    <parameter name="blocks">5</parameter>\n'
        '  </qmc>\n'
        '</simulation>\n'
    )

    info = QmcpackInputInfo(str(filepath))

    assert info.qmc_type == 'dmc'
    assert list(info.qmc_info) == [2, 3, 4, 5, 6]
    assert [qmc.qmc for qmc in info.qmc_info.values()] == [
        'opt', 'vmc', 'opt', 'vmc', 'dmc'
    ]
    assert [qmc.blocks for qmc in info.qmc_info.values()] == [3, 4, 3, 4, 5]


@pytest.mark.parametrize('method', ['linear', 'cslinear', 'linear_batch'])
def test_qmcpack_input_info_optimization_method_aliases(tmp_path, method):
    from ..qmcpack_analyzer_new import QmcpackInputInfo

    filepath = tmp_path / f'{method}.in.xml'
    filepath.write_text(
        '<simulation>\n'
        '  <project id="opt_method" series="0" />\n'
        f'  <qmc method="{method}" />\n'
        '</simulation>\n'
    )

    info = QmcpackInputInfo(str(filepath))

    assert info.qmc_type == 'opt'
    assert info.qmc_info[0].qmc == 'opt'


def test_qmcpack_input_info_default_project(tmp_path):
    from ..qmcpack_analyzer_new import QmcpackInputInfo

    filepath = tmp_path / 'no_project.in.xml'
    filepath.write_text(
        '<simulation>\n'
        '  <qmc method="vmc">\n'
        '    <parameter name="blocks">7</parameter>\n'
        '  </qmc>\n'
        '</simulation>\n'
    )

    info = QmcpackInputInfo(str(filepath))

    assert info.prefix == 'default_project'
    assert info.series_start == 0
    assert info.qmc_type == 'vmc'
    assert info.qmc_info[0].outfiles == (
        'default_project.s000.scalar.dat',
    )


def test_qmcpack_input_info_projectless_input_fixture():
    from ..qmcpack_analyzer_new import QmcpackInputInfo

    info = QmcpackInputInfo(str(INPUT_FILES / 'OH_mixed_pos.in.xml'))

    assert info.prefix == 'default_project'
    assert info.series_start == 0
    assert info.qmc_type is None
    assert len(info.qmc_info) == 0
