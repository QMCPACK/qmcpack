import pytest
from . import NexusTestOrder
pytestmark = pytest.mark.order(NexusTestOrder.QMCA)

from ..generic import generic_settings
generic_settings.raise_error = True

import sys
import os
from pathlib import Path
from . import TEST_DIR
from ..testing import execute,text_eq


QMCA_EXE = TEST_DIR.parent / "bin/qmca"

QA_PATHS = {
    "vmc":       TEST_DIR / "test_qmcpack_analyzer_files/diamond_gamma/vmc",
    "opt":       TEST_DIR / "test_qmcpack_analyzer_files/diamond_gamma/opt",
    "dmc":       TEST_DIR / "test_qmcpack_analyzer_files/diamond_gamma/dmc",
    "vmc_twist": TEST_DIR / "test_qmcpack_analyzer_files/diamond_twist/vmc",
    "multi":     TEST_DIR / "test_qmcpack_analyzer_files/diamond_gamma",
    }


def test_help():

    help_text = 'Usage: qmca'

    command = f'{sys.executable} {QMCA_EXE}'
    out,err,rc = execute(command)
    assert(help_text in out)

    command = f'{sys.executable} {QMCA_EXE} -h'
    out,err,rc = execute(command)
    assert(help_text in out)
#end def test_help



def test_examples():

    example_text = 'QMCA examples'

    command = f'{sys.executable} {QMCA_EXE} -x'
    out,err,rc = execute(command)
    assert(example_text in out)
#end def test_examples



def test_unit_conversion():

    cwd = Path.cwd()
    os.chdir(QA_PATHS["vmc"])

    command = f'{sys.executable} {QMCA_EXE} -e 5 -q e -u eV --fp=16.8f *scalar*'
    out,err,rc = execute(command)

    out_ref = '''
        vmc  series 0  LocalEnergy  =  -284.62368087 +/- 0.10344802
        '''

    assert(text_eq(out,out_ref))

    os.chdir(cwd)
#end def test_unit_conversion



def test_selected_quantities():

    cwd = Path.cwd()
    os.chdir(QA_PATHS["vmc"])

    command = f"{sys.executable} {QMCA_EXE} -e 5 -q 'e k p' --fp=16.8f *scalar*"
    out,err,rc = execute(command)

    out_ref = '''
        vmc  series 0 
          LocalEnergy           =      -10.45972798 +/-       0.00380164 
          Kinetic               =       11.08225203 +/-       0.02281909 
          LocalPotential        =      -21.54198001 +/-       0.02390850 
        '''

    assert(text_eq(out,out_ref))

    os.chdir(cwd)
#end def test_selected_quantities



def test_all_quantities():

    cwd = Path.cwd()
    os.chdir(QA_PATHS["vmc"])

    command = f"{sys.executable} {QMCA_EXE} -e 5 --fp=16.8f *scalar*"
    out,err,rc = execute(command)

    out_ref = '''
        vmc  series 0 
          LocalEnergy           =      -10.45972798 +/-       0.00380164 
          Variance              =        0.39708591 +/-       0.00971200 
          Kinetic               =       11.08225203 +/-       0.02281909 
          LocalPotential        =      -21.54198001 +/-       0.02390850 
          ElecElec              =       -2.73960796 +/-       0.00627045 
          LocalECP              =       -6.55900348 +/-       0.02845221 
          NonLocalECP           =        0.53229890 +/-       0.00924772 
          IonIon                =      -12.77566747 +/-       0.00000000 
          LocalEnergy_sq        =      109.80363374 +/-       0.07821849 
          MPC                   =       -2.47615080 +/-       0.00644218 
          BlockWeight           =      600.00000000 +/-       0.00000000 
          BlockCPU              =        0.03578981 +/-       0.00022222 
          AcceptRatio           =        0.76940789 +/-       0.00038843 
          Efficiency            =   122102.67468280 +/-       0.00000000 
          TotalTime             =        3.40003187 +/-       0.00000000 
          TotalSamples          =    57000.00000000 +/-       0.00000000 
          -------------------------------------------------------------- 
          CorrectedEnergy       =      -10.19627082 +/-       0.00469500 
        '''

    assert(text_eq(out,out_ref))

    os.chdir(cwd)
#end def test_all_quantities



def test_energy_variance():

    cwd = Path.cwd()
    os.chdir(QA_PATHS["opt"])

    command = f"{sys.executable} {QMCA_EXE} -e 5 -q ev --fp=16.8f *scalar*"
    out,err,rc = execute(command)

    out_ref = '''
              LocalEnergy                 Variance                   ratio 
opt series 0 -10.44922550 +/- 0.00475437  0.60256216 +/- 0.01476264  0.0577 
opt series 1 -10.45426389 +/- 0.00320561  0.40278743 +/- 0.01716415  0.0385 
opt series 2 -10.45991696 +/- 0.00271802  0.39865602 +/- 0.00934316  0.0381 
opt series 3 -10.45830307 +/- 0.00298106  0.38110459 +/- 0.00529809  0.0364 
opt series 4 -10.46298481 +/- 0.00561322  0.38927957 +/- 0.01204068  0.0372 
opt series 5 -10.46086055 +/- 0.00375811  0.39354343 +/- 0.00913372  0.0376 
        '''

    assert(text_eq(out,out_ref))

    os.chdir(cwd)
#end def test_energy_variance



def test_multiple_equilibration():

    cwd = Path.cwd()
    os.chdir(QA_PATHS["dmc"])

    command = f"{sys.executable} {QMCA_EXE} -e '5 10 15 20' -q ev --fp=16.8f *scalar*"
    out,err,rc = execute(command)

    out_ref = '''
              LocalEnergy                 Variance                   ratio 
dmc series 0 -10.48910618 +/- 0.00714379  0.45789123 +/- 0.04510618  0.0437 
dmc series 1 -10.53167630 +/- 0.00163992  0.38531028 +/- 0.00162825  0.0366 
dmc series 2 -10.53061971 +/- 0.00123791  0.38172188 +/- 0.00124608  0.0362 
dmc series 3 -10.52807733 +/- 0.00122687  0.38565052 +/- 0.00196074  0.0366 
        '''

    assert(text_eq(out,out_ref))

    os.chdir(cwd)
#end def test_multiple_equilibration



def test_join():

    cwd = Path.cwd()
    os.chdir(QA_PATHS["dmc"])

    command = f"{sys.executable} {QMCA_EXE} -e 5 -j '1 3' -q ev --fp=16.8f *scalar*"
    out,err,rc = execute(command)

    out_ref = '''
               LocalEnergy                 Variance                   ratio 
dmc  series 0 -10.48910618 +/- 0.00714379  0.45789123 +/- 0.04510618  0.0437 
dmc  series 1 -10.53022752 +/- 0.00073527  0.38410495 +/- 0.00082972  0.0365 
        '''

    assert(text_eq(out,out_ref))

    os.chdir(cwd)
#end def test_join



def test_multiple_directories():

    cwd = Path.cwd()
    os.chdir(QA_PATHS["multi"])

    command = f"{sys.executable} {QMCA_EXE} -e 5 -q ev --fp=16.8f */*scalar*"
    out,err,rc = execute(command)

    out_ref = '''
                  LocalEnergy                 Variance                   ratio 
dmc/dmc series 0 -10.48910618 +/- 0.00714379  0.45789123 +/- 0.04510618  0.0437 
dmc/dmc series 1 -10.53165002 +/- 0.00153762  0.38502189 +/- 0.00158532  0.0366 
dmc/dmc series 2 -10.53018917 +/- 0.00106340  0.38161141 +/- 0.00111391  0.0362 
dmc/dmc series 3 -10.52884336 +/- 0.00135513  0.38568155 +/- 0.00151082  0.0366 
 
opt/opt series 0 -10.44922550 +/- 0.00475437  0.60256216 +/- 0.01476264  0.0577 
opt/opt series 1 -10.45426389 +/- 0.00320561  0.40278743 +/- 0.01716415  0.0385 
opt/opt series 2 -10.45991696 +/- 0.00271802  0.39865602 +/- 0.00934316  0.0381 
opt/opt series 3 -10.45830307 +/- 0.00298106  0.38110459 +/- 0.00529809  0.0364 
opt/opt series 4 -10.46298481 +/- 0.00561322  0.38927957 +/- 0.01204068  0.0372 
opt/opt series 5 -10.46086055 +/- 0.00375811  0.39354343 +/- 0.00913372  0.0376 

vmc/vmc series 0 -10.45972798 +/- 0.00380164  0.39708591 +/- 0.00971200  0.0380
        '''

    assert(text_eq(out,out_ref))

    os.chdir(cwd)
#end def test_multiple_directories



def test_twist_average():

    cwd = Path.cwd()
    os.chdir(QA_PATHS["vmc_twist"])

    command = f"{sys.executable} {QMCA_EXE} -a -e 5 -q ev --fp=16.8f *scalar*"
    out,err,rc = execute(command)

    out_ref = '''
              LocalEnergy                 Variance                   ratio 
avg series 0 -11.34367335 +/- 0.00257603  0.57340688 +/- 0.00442552  0.0505 
        '''

    assert(text_eq(out,out_ref))

    os.chdir(cwd)
#end def test_twist_average



def test_weighted_twist_average():

    cwd = Path.cwd()
    os.chdir(QA_PATHS["vmc_twist"])

    command = f"{sys.executable} {QMCA_EXE} -a -w '1 3 3 1' -e 5 -q ev --fp=16.8f *scalar*"
    out,err,rc = execute(command)

    out_ref = '''
              LocalEnergy                 Variance                   ratio 
avg series 0 -11.44375840 +/- 0.00292164  0.44863011 +/- 0.00502859  0.0392 
        '''

    assert(text_eq(out,out_ref))

    os.chdir(cwd)
#end def test_weighted_twist_average



def test_save_plot_pdf():

    try:
        import matplotlib.pyplot  # noqa: F401
    except ImportError:
        pytest.skip('matplotlib not available')

    cwd = Path.cwd()
    os.chdir(QA_PATHS["vmc"])

    prefix = 'testplot_pdf'
    plot_file = f'{prefix}_vmc_LocalEnergy_trace.pdf'
    if os.path.exists(plot_file):
        os.remove(plot_file)

    command = (f'{sys.executable} {QMCA_EXE} -t -q e -e 5 '
               f'--image-prefix {prefix} --image-format pdf --nowarn *scalar*')
    out,err,rc = execute(command)
    assert rc==0
    assert os.path.exists(plot_file)

    os.remove(plot_file)
    os.chdir(cwd)
#end def test_save_plot_pdf



def test_save_plot_png_default():

    try:
        import matplotlib.pyplot  # noqa: F401
    except ImportError:
        pytest.skip('matplotlib not available')

    cwd = Path.cwd()
    os.chdir(QA_PATHS["vmc"])

    prefix = 'testplot_png'
    plot_file = f'{prefix}_vmc_LocalEnergy_trace.png'
    if os.path.exists(plot_file):
        os.remove(plot_file)

    command = (f'{sys.executable} {QMCA_EXE} -t -q e -e 5 '
               f'--image-prefix {prefix} --image --nowarn *scalar*')
    out,err,rc = execute(command)
    assert rc==0
    assert os.path.exists(plot_file)

    os.remove(plot_file)
    os.chdir(cwd)
#end def test_save_plot_png_default

