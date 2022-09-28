
import sys
import testing
from testing import execute,text_eq


test_info = dict()
directories = dict()


def get_test_info():
    if len(test_info)==0:
        import os

        tpath = testing.setup_unit_test_output_directory('qmca','test_qmca')
        
        exe = testing.executable_path('qmca')

        dirs = ['diamond_gamma','diamond_twist']

        for dir in dirs:
            qa_files_path = testing.unit_test_file_path('qmcpack_analyzer',dir)
            command = 'rsync -a {} {}'.format(qa_files_path,tpath)
            out,err,rc = execute(command)
            assert(rc==0)
            data_path = os.path.join(tpath,dir)
            assert(os.path.exists(data_path))
        #end for

        paths = dict(
            vmc       = 'diamond_gamma/vmc',
            opt       = 'diamond_gamma/opt',
            dmc       = 'diamond_gamma/dmc',
            vmc_twist = 'diamond_twist/vmc',
            multi     = 'diamond_gamma',
            )

        test_info['tpath'] = tpath
        test_info['exe']   = exe
        for name,path in paths.items():
            directories[name] = os.path.join(tpath,path)
        #end for
    #end if
#end def get_test_info


def get_exe():
    get_test_info()
    return test_info['exe']
#end def get_exe    


def enter(path_name):
    import os
    assert(path_name in directories)
    path = directories[path_name]
    assert(os.path.exists(path))
    assert('cwd' not in directories)
    directories['cwd'] = os.getcwd()
    assert('cwd' in directories)
    os.chdir(path)
#end def enter


def leave():
    import os
    assert('cwd' in directories)
    cwd = directories.pop('cwd')
    assert(not 'cwd' in directories)
    os.chdir(cwd)
#end def leave



def test_help():
    exe = get_exe()

    help_text = 'Usage: qmca'

    command = sys.executable+' {}'.format(exe)
    out,err,rc = execute(command)
    assert(help_text in out)

    command = sys.executable+' {} -h'.format(exe)
    out,err,rc = execute(command)
    assert(help_text in out)
#end def test_help



def test_examples():
    exe = get_exe()

    example_text = 'QMCA examples'

    command = sys.executable+' {} -x'.format(exe)
    out,err,rc = execute(command)
    assert(example_text in out)
#end def test_examples



def test_unit_conversion():
    exe = get_exe()

    enter('vmc')

    command = sys.executable+' {} -e 5 -q e -u eV --fp=16.8f *scalar*'.format(exe)
    out,err,rc = execute(command)

    out_ref = '''
        vmc  series 0  LocalEnergy  =  -284.62368087 +/- 0.10344802
        '''

    assert(text_eq(out,out_ref))

    leave()
#end def test_unit_conversion



def test_selected_quantities():
    exe = get_exe()

    enter('vmc')

    command = sys.executable+" {} -e 5 -q 'e k p' --fp=16.8f *scalar*".format(exe)
    out,err,rc = execute(command)

    out_ref = '''
        vmc  series 0 
          LocalEnergy           =      -10.45972798 +/-       0.00380164 
          Kinetic               =       11.08225203 +/-       0.02281909 
          LocalPotential        =      -21.54198001 +/-       0.02390850 
        '''

    assert(text_eq(out,out_ref))

    leave()
#end def test_selected_quantities



def test_all_quantities():
    exe = get_exe()

    enter('vmc')

    command = sys.executable+" {} -e 5 --fp=16.8f *scalar*".format(exe)
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

    leave()
#end def test_all_quantities



def test_energy_variance():
    exe = get_exe()

    enter('opt')

    command = sys.executable+" {} -e 5 -q ev --fp=16.8f *scalar*".format(exe)
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

    leave()
#end def test_energy_variance



def test_multiple_equilibration():
    exe = get_exe()

    enter('dmc')

    command = sys.executable+" {} -e '5 10 15 20' -q ev --fp=16.8f *scalar*".format(exe)
    out,err,rc = execute(command)

    out_ref = '''
              LocalEnergy                 Variance                   ratio 
dmc series 0 -10.48910618 +/- 0.00714379  0.45789123 +/- 0.04510618  0.0437 
dmc series 1 -10.53167630 +/- 0.00163992  0.38531028 +/- 0.00162825  0.0366 
dmc series 2 -10.53061971 +/- 0.00123791  0.38172188 +/- 0.00124608  0.0362 
dmc series 3 -10.52807733 +/- 0.00122687  0.38565052 +/- 0.00196074  0.0366 
        '''

    assert(text_eq(out,out_ref))

    leave()
#end def test_multiple_equilibration



def test_join():
    exe = get_exe()

    enter('dmc')

    command = sys.executable+" {} -e 5 -j '1 3' -q ev --fp=16.8f *scalar*".format(exe)
    out,err,rc = execute(command)

    out_ref = '''
               LocalEnergy                 Variance                   ratio 
dmc  series 0 -10.48910618 +/- 0.00714379  0.45789123 +/- 0.04510618  0.0437 
dmc  series 1 -10.53022752 +/- 0.00073527  0.38410495 +/- 0.00082972  0.0365 
        '''

    assert(text_eq(out,out_ref))

    leave()
#end def test_join



def test_multiple_directories():
    exe = get_exe()

    enter('multi')

    command = sys.executable+" {} -e 5 -q ev --fp=16.8f */*scalar*".format(exe)
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

    leave()
#end def test_multiple_directories



def test_twist_average():
    exe = get_exe()

    enter('vmc_twist')

    command = sys.executable+" {} -a -e 5 -q ev --fp=16.8f *scalar*".format(exe)
    out,err,rc = execute(command)

    out_ref = '''
              LocalEnergy                 Variance                   ratio 
avg series 0 -11.34367335 +/- 0.00257603  0.57340688 +/- 0.00442552  0.0505 
        '''

    assert(text_eq(out,out_ref))

    leave()
#end def test_twist_average



def test_weighted_twist_average():
    exe = get_exe()

    enter('vmc_twist')

    command = sys.executable+" {} -a -w '1 3 3 1' -e 5 -q ev --fp=16.8f *scalar*".format(exe)
    out,err,rc = execute(command)

    out_ref = '''
              LocalEnergy                 Variance                   ratio 
avg series 0 -11.44375840 +/- 0.00292164  0.44863011 +/- 0.00502859  0.0392 
        '''

    assert(text_eq(out,out_ref))

    leave()
#end def test_weighted_twist_average

