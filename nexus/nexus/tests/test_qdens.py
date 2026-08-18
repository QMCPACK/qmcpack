import shutil

import numpy as np
import pytest

h5py = pytest.importorskip('h5py')

from . import NexusTestOrder
pytestmark = pytest.mark.order(NexusTestOrder.QDENS)



from . import TEST_DIR
from ..testing import execute, text_eq


def test_density(tmp_path):

    exe = TEST_DIR.parent / "bin/qdens"

    qa_files_path = TEST_DIR / "test_qmcpack_analyzer_files/diamond_gamma"
    shutil.copytree(qa_files_path, tmp_path, dirs_exist_ok=True)

    dmc_path = tmp_path / 'dmc'

    files_bef = (
        dmc_path / "dmc.s000.stat.h5",
        dmc_path / "dmc.s000.scalar.dat",
        dmc_path / "dmc.s001.stat.h5",
        dmc_path / "dmc.s001.scalar.dat",
        dmc_path / "dmc.s002.stat.h5",
        dmc_path / "dmc.s002.scalar.dat",
        dmc_path / "dmc.s003.stat.h5",
        dmc_path / "dmc.s003.scalar.dat",
        dmc_path / "dmc.in.xml",
        dmc_path / "dmc.out",
        dmc_path / "dmc.err",
        )

    assert(set(dmc_path.iterdir()) == set(files_bef))

    command = f'{exe} -v -e 4 -f xsf -i {dmc_path}/dmc.in.xml {dmc_path}/*stat.h5'

    out,err,rc = execute(command)

    files_aft = (
        tmp_path / "dmc/dmc.err",
        tmp_path / "dmc/dmc.in.xml",
        tmp_path / "dmc/dmc.out",
        tmp_path / "dmc/dmc.s000.scalar.dat",
        tmp_path / "dmc/dmc.s000.SpinDensity_d+err.xsf",
        tmp_path / "dmc/dmc.s000.SpinDensity_d-err.xsf",
        tmp_path / "dmc/dmc.s000.SpinDensity_d.xsf",
        tmp_path / "dmc/dmc.s000.SpinDensity_u+d+err.xsf",
        tmp_path / "dmc/dmc.s000.SpinDensity_u+d-err.xsf",
        tmp_path / "dmc/dmc.s000.SpinDensity_u-d+err.xsf",
        tmp_path / "dmc/dmc.s000.SpinDensity_u-d-err.xsf",
        tmp_path / "dmc/dmc.s000.SpinDensity_u+d.xsf",
        tmp_path / "dmc/dmc.s000.SpinDensity_u-d.xsf",
        tmp_path / "dmc/dmc.s000.SpinDensity_u+err.xsf",
        tmp_path / "dmc/dmc.s000.SpinDensity_u-err.xsf",
        tmp_path / "dmc/dmc.s000.SpinDensity_u.xsf",
        tmp_path / "dmc/dmc.s000.stat.h5",
        tmp_path / "dmc/dmc.s001.scalar.dat",
        tmp_path / "dmc/dmc.s001.SpinDensity_d+err.xsf",
        tmp_path / "dmc/dmc.s001.SpinDensity_d-err.xsf",
        tmp_path / "dmc/dmc.s001.SpinDensity_d.xsf",
        tmp_path / "dmc/dmc.s001.SpinDensity_u+d+err.xsf",
        tmp_path / "dmc/dmc.s001.SpinDensity_u+d-err.xsf",
        tmp_path / "dmc/dmc.s001.SpinDensity_u-d+err.xsf",
        tmp_path / "dmc/dmc.s001.SpinDensity_u-d-err.xsf",
        tmp_path / "dmc/dmc.s001.SpinDensity_u+d.xsf",
        tmp_path / "dmc/dmc.s001.SpinDensity_u-d.xsf",
        tmp_path / "dmc/dmc.s001.SpinDensity_u+err.xsf",
        tmp_path / "dmc/dmc.s001.SpinDensity_u-err.xsf",
        tmp_path / "dmc/dmc.s001.SpinDensity_u.xsf",
        tmp_path / "dmc/dmc.s001.stat.h5",
        tmp_path / "dmc/dmc.s002.scalar.dat",
        tmp_path / "dmc/dmc.s002.SpinDensity_d+err.xsf",
        tmp_path / "dmc/dmc.s002.SpinDensity_d-err.xsf",
        tmp_path / "dmc/dmc.s002.SpinDensity_d.xsf",
        tmp_path / "dmc/dmc.s002.SpinDensity_u+d+err.xsf",
        tmp_path / "dmc/dmc.s002.SpinDensity_u+d-err.xsf",
        tmp_path / "dmc/dmc.s002.SpinDensity_u-d+err.xsf",
        tmp_path / "dmc/dmc.s002.SpinDensity_u-d-err.xsf",
        tmp_path / "dmc/dmc.s002.SpinDensity_u+d.xsf",
        tmp_path / "dmc/dmc.s002.SpinDensity_u-d.xsf",
        tmp_path / "dmc/dmc.s002.SpinDensity_u+err.xsf",
        tmp_path / "dmc/dmc.s002.SpinDensity_u-err.xsf",
        tmp_path / "dmc/dmc.s002.SpinDensity_u.xsf",
        tmp_path / "dmc/dmc.s002.stat.h5",
        tmp_path / "dmc/dmc.s003.scalar.dat",
        tmp_path / "dmc/dmc.s003.SpinDensity_d+err.xsf",
        tmp_path / "dmc/dmc.s003.SpinDensity_d-err.xsf",
        tmp_path / "dmc/dmc.s003.SpinDensity_d.xsf",
        tmp_path / "dmc/dmc.s003.SpinDensity_u+d+err.xsf",
        tmp_path / "dmc/dmc.s003.SpinDensity_u+d-err.xsf",
        tmp_path / "dmc/dmc.s003.SpinDensity_u-d+err.xsf",
        tmp_path / "dmc/dmc.s003.SpinDensity_u-d-err.xsf",
        tmp_path / "dmc/dmc.s003.SpinDensity_u+d.xsf",
        tmp_path / "dmc/dmc.s003.SpinDensity_u-d.xsf",
        tmp_path / "dmc/dmc.s003.SpinDensity_u+err.xsf",
        tmp_path / "dmc/dmc.s003.SpinDensity_u-err.xsf",
        tmp_path / "dmc/dmc.s003.SpinDensity_u.xsf",
        tmp_path / "dmc/dmc.s003.stat.h5",
        )

    assert(set(dmc_path.iterdir()) == set(files_aft))

    tot_file = dmc_path / 'dmc.s003.SpinDensity_u+d.xsf'
    pol_file = dmc_path / 'dmc.s003.SpinDensity_u-d.xsf'

    tot = tot_file.read_text()
    pol = pol_file.read_text()

    tot_ref = '''
         CRYSTAL
         PRIMVEC 
             1.78500000   1.78500000   0.00000000
            -0.00000000   1.78500000   1.78500000
             1.78500000  -0.00000000   1.78500000
         PRIMCOORD 
           2 1
             6   0.00000000   0.00000000   0.00000000
             6   0.89250000   0.89250000   0.89250000
         BEGIN_BLOCK_DATAGRID_3D
           density
           BEGIN_DATAGRID_3D_density
             4 4 4
             0.59500000   0.59500000   0.59500000
             1.78500000   1.78500000   0.00000000
            -0.00000000   1.78500000   1.78500000
             1.78500000  -0.00000000   1.78500000
               0.73126076   0.62407496   0.51676366   0.73126076
               0.62575089   0.19225114   0.18686389   0.62575089
               0.51847569   0.18457799   0.42203355   0.51847569
               0.73126076   0.62407496   0.51676366   0.73126076
               0.62659840   0.19325900   0.18422995   0.62659840
               0.19219866   0.04873728   0.13184395   0.19219866
               0.18474638   0.13013188   0.10227670   0.18474638
               0.62659840   0.19325900   0.18422995   0.62659840
               0.51793019   0.18615766   0.41806405   0.51793019
               0.18425005   0.13092538   0.10088238   0.18425005
               0.41967003   0.10133434   0.14471118   0.41967003
               0.51793019   0.18615766   0.41806405   0.51793019
               0.73126076   0.62407496   0.51676366   0.73126076
               0.62575089   0.19225114   0.18686389   0.62575089
               0.51847569   0.18457799   0.42203355   0.51847569
               0.73126076   0.62407496   0.51676366   0.73126076
           END_DATAGRID_3D_density
         END_BLOCK_DATAGRID_3D
        '''

    pol_ref = '''
         CRYSTAL
         PRIMVEC 
             1.78500000   1.78500000   0.00000000
            -0.00000000   1.78500000   1.78500000
             1.78500000  -0.00000000   1.78500000
         PRIMCOORD 
           2 1
             6   0.00000000   0.00000000   0.00000000
             6   0.89250000   0.89250000   0.89250000
         BEGIN_BLOCK_DATAGRID_3D
           density
           BEGIN_DATAGRID_3D_density
             4 4 4
             0.59500000   0.59500000   0.59500000
             1.78500000   1.78500000   0.00000000
            -0.00000000   1.78500000   1.78500000
             1.78500000  -0.00000000   1.78500000
               0.00106753   0.00015792  -0.00122859   0.00106753
              -0.00003402   0.00018762  -0.00051347  -0.00003402
               0.00154254   0.00067654   0.00073434   0.00154254
               0.00106753   0.00015792  -0.00122859   0.00106753
               0.00263956   0.00079744  -0.00118289   0.00263956
              -0.00039348  -0.00026396  -0.00069392  -0.00039348
               0.00087929   0.00000719   0.00113934   0.00087929
               0.00263956   0.00079744  -0.00118289   0.00263956
              -0.00013655  -0.00041508  -0.00235212  -0.00013655
               0.00003805  -0.00025962  -0.00133495   0.00003805
               0.00040692   0.00051699  -0.00198263   0.00040692
              -0.00013655  -0.00041508  -0.00235212  -0.00013655
               0.00106753   0.00015792  -0.00122859   0.00106753
              -0.00003402   0.00018762  -0.00051347  -0.00003402
               0.00154254   0.00067654   0.00073434   0.00154254
               0.00106753   0.00015792  -0.00122859   0.00106753
           END_DATAGRID_3D_density
         END_BLOCK_DATAGRID_3D
        '''

    assert(text_eq(tot,tot_ref,atol=1e-7))
    assert(text_eq(pol,pol_ref,atol=1e-7))
#end def test_density


def test_qmc_spindensity_metadata_uses_cell_and_corner(tmp_path):
    exe = TEST_DIR.parent / 'bin/qdens'

    input_file = tmp_path / 'he.xml'
    input_file.write_text('''\
<?xml version="1.0"?>
<simulation>
  <project id="case" series="0"/>
  <qmcsystem>
    <simulationcell>
      <parameter name="lattice" units="bohr">
        5 0 0
        0 5 0
        0 0 5
      </parameter>
      <parameter name="bconds">n n n</parameter>
    </simulationcell>
    <particleset name="ion0" size="1">
      <group name="He"><parameter name="charge">2</parameter></group>
      <attrib name="position" datatype="posArray">0 0 0</attrib>
    </particleset>
    <particleset name="e">
      <group name="u" size="1"><parameter name="charge">-1</parameter></group>
      <group name="d" size="1"><parameter name="charge">-1</parameter></group>
    </particleset>
  </qmcsystem>
  <qmc method="vmc_batch" move="pbyp">
    <estimators>
      <estimator name="spindensity" type="spindensity">
        <parameter name="grid">2 2 2</parameter>
        <parameter name="corner">0 0 0</parameter>
        <parameter name="cell">
          2 0 0
          0 3 0
          0 0 4
        </parameter>
      </estimator>
    </estimators>
  </qmc>
</simulation>
''')

    stat_file = tmp_path / 'case.s000.stat.h5'
    with h5py.File(stat_file,'w') as h5:
        g = h5.create_group('spindensity')
        u = g.create_group('u')
        d = g.create_group('d')
        u.create_dataset('value',data=np.ones((2,8)))
        d.create_dataset('value',data=2*np.ones((2,8)))

    command = f'{exe} -f xsf -i {input_file} {stat_file}'
    execute(command)

    xsf = (tmp_path / 'case.s000.spindensity_u+d.xsf').read_text()
    assert '3 3 3' in xsf
    assert '1.05835442E+00 0.00000000E+00 0.00000000E+00' in xsf
    assert '0.00000000E+00 1.58753163E+00 0.00000000E+00' in xsf
    assert '0.00000000E+00 0.00000000E+00 2.11670883E+00' in xsf

    with pytest.raises(
            AssertionError,
            match='cannot write CHGCAR with a density cell different from the simulation cell'):
        execute(f'{exe} -f chgcar -i {input_file} {stat_file}')
    #end with
#end def test_qmc_spindensity_metadata_uses_cell_and_corner


def _write_series_spindensity_input(filepath):
    filepath.write_text('''\
<?xml version="1.0"?>
<simulation>
  <project id="series_spindensity" series="5"/>
  <estimators>
    <estimator name="GlobalSpinDensity" type="spindensity">
      <parameter name="grid">2 2 2</parameter>
    </estimator>
  </estimators>
  <qmcsystem>
    <simulationcell>
      <parameter name="lattice" units="bohr">
        5 0 0
        0 5 0
        0 0 5
      </parameter>
      <parameter name="bconds">p p p</parameter>
    </simulationcell>
    <particleset name="ion0" size="1">
      <group name="He"><parameter name="charge">2</parameter></group>
      <attrib name="position" datatype="posArray">0 0 0</attrib>
    </particleset>
    <particleset name="e">
      <group name="u" size="1"><parameter name="charge">-1</parameter></group>
      <group name="d" size="1"><parameter name="charge">-1</parameter></group>
    </particleset>
  </qmcsystem>
  <qmc method="vmc_batch" move="pbyp">
    <estimators>
      <estimator name="SpinDensity" type="spindensity">
        <parameter name="grid">2 2 2</parameter>
        <parameter name="corner">0 0 0</parameter>
        <parameter name="cell">2 0 0  0 3 0  0 0 4</parameter>
      </estimator>
    </estimators>
  </qmc>
  <qmc method="dmc_batch" move="pbyp">
    <estimators>
      <estimator name="SpinDensity" type="spindensity">
        <parameter name="grid">2 3 4</parameter>
        <parameter name="corner">1 2 3</parameter>
        <parameter name="cell">4 0 0  0 5 0  0 0 6</parameter>
      </estimator>
    </estimators>
  </qmc>
</simulation>
''')
#end def _write_series_spindensity_input


def _write_spin_density_stat(filepath,spin_density_cells,*,global_cells=8,name='SpinDensity'):
    estimators = [(name,spin_density_cells)]
    if global_cells is not None:
        estimators.insert(0,('GlobalSpinDensity',global_cells))
    #end if
    with h5py.File(filepath,'w') as h5:
        for estimator,cells in estimators:
            group = h5.create_group(estimator)
            group.create_group('u').create_dataset('value',data=np.ones((2,cells)))
            group.create_group('d').create_dataset('value',data=2*np.ones((2,cells)))
#end def _write_spin_density_stat


def test_bounded_density_metadata_uses_fractional_delta(tmp_path):
    exe = TEST_DIR.parent / 'bin/qdens'
    input_file = tmp_path / 'density.xml'
    input_file.write_text('''\
<?xml version="1.0"?>
<simulation>
  <project id="density" series="0"/>
  <qmcsystem>
    <simulationcell>
      <parameter name="lattice" units="bohr">
        6 0 0
        0 6 0
        0 0 6
      </parameter>
      <parameter name="bconds">n n n</parameter>
    </simulationcell>
    <particleset name="ion0" size="1">
      <group name="He"><parameter name="charge">2</parameter></group>
      <attrib name="position" datatype="posArray">0 0 0</attrib>
    </particleset>
    <particleset name="e">
      <group name="u" size="1"><parameter name="charge">-1</parameter></group>
      <group name="d" size="1"><parameter name="charge">-1</parameter></group>
    </particleset>
    <hamiltonian name="h0" type="generic" target="e">
      <estimator name="Density" type="density" delta="0.5 0.25 0.25"
                 x_min="1" x_max="3" y_min="2" y_max="5" z_min="0" z_max="4"/>
    </hamiltonian>
  </qmcsystem>
  <qmc method="vmc"><parameter name="blocks">2</parameter></qmc>
</simulation>
''')
    stat_file = tmp_path / 'density.s000.stat.h5'
    with h5py.File(stat_file,'w') as h5:
        h5.create_group('Density').create_dataset(
            'value',data=np.ones((2,2,4,4),dtype=float))

    execute(f'{exe} -f xsf -i {input_file} {stat_file}')
    xsf = (tmp_path / 'density.s000.Density_q.xsf').read_text()
    assert '3 5 5' in xsf
    assert '7.93765813E-01 1.25679587E+00 2.64588604E-01' in xsf
    assert '1.05835442E+00 0.00000000E+00 0.00000000E+00' in xsf
    assert '0.00000000E+00 1.58753163E+00 0.00000000E+00' in xsf
    assert '0.00000000E+00 0.00000000E+00 2.11670883E+00' in xsf
#end def test_bounded_density_metadata_uses_fractional_delta


@pytest.mark.parametrize(
    'lattice,bconds,delta,message',
    (
        ((6,0,0,0,6,0,0,0,6),'p p p','0.5 0.25 0.25',
         'open boundary conditions'),
        ((6,1,0,0,6,0,0,0,6),'n n n','0.5 0.25 0.25',
         'orthorhombic simulation cell'),
        ((6,0,0,0,6,0,0,0,6),'n n n','0 0.25 0.25',
         'three positive delta values'),
        ),
    )
def test_bounded_density_metadata_validates_cell_and_delta(
        tmp_path,lattice,bconds,delta,message):
    lattice = '\n        '.join(
        ' '.join(str(value) for value in axis) for axis in np.array(lattice).reshape(3,3))
    exe = TEST_DIR.parent / 'bin/qdens'
    input_file = tmp_path / 'density.xml'
    input_file.write_text(f'''\
<?xml version="1.0"?>
<simulation>
  <project id="density" series="0"/>
  <qmcsystem>
    <simulationcell>
      <parameter name="lattice" units="bohr">
        {lattice}
      </parameter>
      <parameter name="bconds">{bconds}</parameter>
    </simulationcell>
    <particleset name="ion0" size="1">
      <group name="He"><parameter name="charge">2</parameter></group>
      <attrib name="position" datatype="posArray">0 0 0</attrib>
    </particleset>
    <particleset name="e">
      <group name="u" size="1"><parameter name="charge">-1</parameter></group>
      <group name="d" size="1"><parameter name="charge">-1</parameter></group>
    </particleset>
    <hamiltonian name="h0" type="generic" target="e">
      <estimator name="Density" type="density" delta="{delta}"
                 x_min="1" x_max="3" y_min="2" y_max="5" z_min="0" z_max="4"/>
    </hamiltonian>
  </qmcsystem>
  <qmc method="vmc"><parameter name="blocks">2</parameter></qmc>
</simulation>
''')
    stat_file = tmp_path / 'density.s000.stat.h5'
    with h5py.File(stat_file,'w') as h5:
        h5.create_group('Density').create_dataset(
            'value',data=np.ones((2,2,4,4),dtype=float))

    with pytest.raises(AssertionError, match=message):
        execute(f'{exe} -f xsf -i {input_file} {stat_file}')
    #end with
#end def test_bounded_density_metadata_validates_cell_and_delta


def test_spindensity_metadata_is_selected_per_series(tmp_path):
    exe = TEST_DIR.parent / 'bin/qdens'
    infile = tmp_path / 'series_spindensity.xml'
    _write_series_spindensity_input(infile)
    stat5 = tmp_path / 'series_spindensity.s005.stat.h5'
    stat6 = tmp_path / 'series_spindensity.s006.stat.h5'
    _write_spin_density_stat(stat5,8)
    _write_spin_density_stat(stat6,24)

    execute(f'{exe} -f xsf -i {infile} {stat5} {stat6}')
    assert (tmp_path / 'series_spindensity.s005.GlobalSpinDensity_u+d.xsf').exists()
    assert (tmp_path / 'series_spindensity.s006.GlobalSpinDensity_u+d.xsf').exists()

    first = (tmp_path / 'series_spindensity.s005.SpinDensity_u+d.xsf').read_text()
    second = (tmp_path / 'series_spindensity.s006.SpinDensity_u+d.xsf').read_text()
    assert '3 3 3' in first
    assert '1.05835442E+00 0.00000000E+00 0.00000000E+00' in first
    assert '0.00000000E+00 1.58753163E+00 0.00000000E+00' in first
    assert '2.64588604E-01 3.96882906E-01 5.29177209E-01' in first
    assert '3 4 5' in second
    assert '2.11670883E+00 0.00000000E+00 0.00000000E+00' in second
    assert '0.00000000E+00 2.64588604E+00 0.00000000E+00' in second
    assert '1.05835442E+00 1.49933542E+00 1.98441453E+00' in second
#end def test_spindensity_metadata_is_selected_per_series


def test_spindensity_input_metadata_overrides_grid_fallback(tmp_path):
    exe = TEST_DIR.parent / 'bin/qdens'
    infile = tmp_path / 'series_spindensity.xml'
    _write_series_spindensity_input(infile)
    stat = tmp_path / 'series_spindensity.s005.stat.h5'
    _write_spin_density_stat(stat,8)

    out,err,_ = execute(f'{exe} -f xsf -g "1 1 8" -i {infile} {stat}')
    assert 'Ignoring --grid' in out + err
    assert '3 3 3' in (tmp_path / 'series_spindensity.s005.SpinDensity_u+d.xsf').read_text()
#end def test_spindensity_input_metadata_overrides_grid_fallback


def test_spindensity_input_hdf5_grid_mismatch_fails_clearly(tmp_path):
    exe = TEST_DIR.parent / 'bin/qdens'
    infile = tmp_path / 'series_spindensity.xml'
    _write_series_spindensity_input(infile)
    stat = tmp_path / 'series_spindensity.s005.stat.h5'
    _write_spin_density_stat(stat,7)

    with pytest.raises(
            AssertionError,
            match='density grid does not match number of HDF5 data cells') as error:
        execute(f'{exe} -f xsf -i {infile} {stat}')
    #end with
    message = str(error.value)
    assert 'series: 5' in message
    assert 'estimator: SpinDensity' in message
    assert 'HDF5 data cells: 7' in message
    assert str(stat) in message
#end def test_spindensity_input_hdf5_grid_mismatch_fails_clearly


def test_spindensity_series_mismatch_requires_unambiguous_fallback(tmp_path):
    exe = TEST_DIR.parent / 'bin/qdens'
    infile = tmp_path / 'series_spindensity.xml'
    _write_series_spindensity_input(infile)
    stat = tmp_path / 'series_spindensity.s007.stat.h5'
    _write_spin_density_stat(stat,8,global_cells=None)

    with pytest.raises(
            AssertionError,
            match='could not find input metadata for density estimator "SpinDensity" in series 7'):
        execute(f'{exe} -f xsf -i {infile} {stat}')
    #end with

    with pytest.raises(
            AssertionError,
            match='density grid does not match number of HDF5 data cells'):
        execute(f'{exe} -f xsf -g "3 3 3" -i {infile} {stat}')
    #end with
#end def test_spindensity_series_mismatch_requires_unambiguous_fallback
