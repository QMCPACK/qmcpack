
import testing
from testing import divert_nexus,restore_nexus
from testing import divert_nexus_log,restore_nexus_log
from testing import value_eq,object_eq


associated_files = dict()

def get_files():
    return testing.collect_unit_test_file_paths('gamess_input',associated_files)
#end def get_files



def make_serial_reference(gi):
    s = gi.serial()
    ref = '    ref = {\n'
    for k in sorted(s.keys()):
        v = s[k]
        if isinstance(v,str):
            v = "'"+v+"'"
        #end if
        ref +="        '{}' : {},\n".format(k,v)
    #end for
    ref += '        }\n'
    return ref
#end def make_serial_reference


serial_references = dict()


h2o_basis_text = '''A molecule.
Cnv 2

O 8       0.00000000       0.00000000       0.00000000
s 9 1.00
1 0.125346     0.055741
2 0.268022     0.304848
3 0.573098     0.453752
4 1.225429     0.295926
5 2.620277     0.019567
6 5.602818     -0.128627
7 11.980245     0.012024
8 25.616801     0.000407
9 54.775216     -0.000076
s 1 1.00
1 0.160664     1.000000
s 1 1.00
1 0.384526     1.000000
s 1 1.00
1 0.935157     1.000000
s 1 1.00
1 1.937532     1.000000
p 9 1.00
1 0.083598     0.044958
2 0.167017     0.150175
3 0.333673     0.255999
4 0.666627     0.281879
5 1.331816     0.242835
6 2.660761     0.161134
7 5.315785     0.082308
8 10.620108     0.039899
9 21.217318     0.004679
p 1 1.00
1 0.130580     1.000000
p 1 1.00
1 0.372674     1.000000
p 1 1.00
1 1.178227     1.000000
p 1 1.00
1 1.589967     1.000000
d 1 1.00
1 0.401152     1.000000
d 1 1.00
1 1.174596     1.000000
d 1 1.00
1 2.823972     1.000000
d 1 1.00
1 4.292433     1.000000
f 1 1.00
1 0.708666     1.000000
f 1 1.00
1 2.006788     1.000000
f 1 1.00
1 3.223721     1.000000
g 1 1.00
1 1.207657     1.000000
g 1 1.00
1 3.584495     1.000000

H 1       0.00000000       0.75716000       0.58626000
s 6 1.00
1 6.359201     -0.004943
2 3.546637     0.049579
3 1.493442     0.037176
4 0.551948     0.287908
5 0.207203     0.009543
6 0.079234     0.770084
s 6 1.00
1 6.359201     -0.016672
2 3.546637     -0.005774
3 1.493442     -0.227982
4 0.551948     -0.285652
5 0.207203     -1.071579
6 0.079234     1.423767
s 6 1.00
1 6.359201     -0.018886
2 3.546637     -0.058854
3 1.493442     -0.556988
4 0.551948     -1.084022
5 0.207203     2.721525
6 0.079234     -1.458091
s 6 1.00
1 6.359201     -0.081999
2 3.546637     -0.069460
3 1.493442     -2.075555
4 0.551948     4.125180
5 0.207203     -3.407591
6 0.079234     1.264873
s 6 1.00
1 6.359201     0.624292
2 3.546637     -3.954593
3 1.493442     5.614244
4 0.551948     -4.299981
5 0.207203     2.379542
6 0.079234     -0.749499
s 1 1.00
1 0.102700     1.000000
p 1 1.00
1 4.516000     1.000000
p 1 1.00
1 1.712000     1.000000
p 1 1.00
1 0.649000     1.000000
p 1 1.00
1 0.246000     1.000000
d 1 1.00
1 2.950000     1.000000
d 1 1.00
1 1.206000     1.000000
d 1 1.00
1 0.493000     1.000000
f 1 1.00
1 1.397000     1.000000
f 1 1.00
1 2.506000     1.000000
f 1 1.00
1 0.875000     1.000000
g 1 1.00
1 2.358000     1.000000

'''

h2o_ecp_text = '''O-QMC GEN 2 1
3
6.00000000 1 9.29793903
55.78763416 3 8.86492204
-38.81978498 2 8.62925665
1
38.41914135 2 8.71924452
H-QMC GEN 0 0
3
1.00000000 1 4.47692410
4.47692410 3 2.97636451
-4.32112340 2 3.01841596
H-QMC
'''



def generate_serial_references():

    basis_text = h2o_basis_text
    ecp_text   = h2o_ecp_text

    serial_references['rhf.inp'] = {
        'contrl/coord'  : 'unique',
        'contrl/ecp'    : 'read',
        'contrl/exetyp' : 'run',
        'contrl/icharg' : 0,
        'contrl/ispher' : 1,
        'contrl/maxit'  : 200,
        'contrl/mult'   : 1,
        'contrl/runtyp' : 'energy',
        'contrl/scftyp' : 'rohf',
        'data/text'     : basis_text,
        'ecp/text'      : ecp_text,
        'guess/guess'   : 'huckel',
        'scf/dirscf'    : True,
        'system/memory' : 150000000,
        }

    serial_references['cisd.inp'] = {
        'cidrt/group'   : 'c2v',
        'cidrt/iexcit'  : 2,
        'cidrt/istsym'  : 1,
        'cidrt/mxnint'  : 500000,
        'cidrt/nalp'    : 0,
        'cidrt/ndoc'    : 4,
        'cidrt/nfzc'    : 0,
        'cidrt/nprt'    : 2,
        'cidrt/nval'    : 60,
        'contrl/cityp'  : 'guga',
        'contrl/coord'  : 'unique',
        'contrl/ecp'    : 'read',
        'contrl/exetyp' : 'run',
        'contrl/icharg' : 0,
        'contrl/ispher' : 1,
        'contrl/maxit'  : 200,
        'contrl/mult'   : 1,
        'contrl/runtyp' : 'energy',
        'contrl/scftyp' : 'none',
        'data/text'     : basis_text,
        'ecp/text'      : ecp_text,
        'gugdia/cvgtol' : 1e-05,
        'gugdia/itermx' : 1000,
        'gugdia/prttol' : 0.001,
        'scf/dirscf'    : True,
        'system/memory' : 150000000,
        }

    serial_references['cas.inp'] = {
        'contrl/coord'  : 'unique',
        'contrl/ecp'    : 'read',
        'contrl/exetyp' : 'run',
        'contrl/icharg' : 0,
        'contrl/ispher' : 1,
        'contrl/maxit'  : 200,
        'contrl/mult'   : 1,
        'contrl/runtyp' : 'energy',
        'contrl/scftyp' : 'mcscf',
        'data/text'     : basis_text,
        'drt/fors'      : True,
        'drt/group'     : 'c2v',
        'drt/istsym'    : 1,
        'drt/mxnint'    : 500000,
        'drt/nalp'      : 0,
        'drt/ndoc'      : 4,
        'drt/nmcc'      : 0,
        'drt/nval'      : 4,
        'ecp/text'      : ecp_text,
        'mcscf/acurcy'  : 1e-05,
        'mcscf/cistep'  : 'guga',
        'mcscf/fullnr'  : True,
        'mcscf/maxit'   : 1000,
        'scf/dirscf'    : True,
        'system/memory' : 150000000,
        }

#end def generate_serial_references


def get_serial_references():
    if len(serial_references)==0:
        generate_serial_references()
    #end if
    return serial_references
#end def get_serial_references


def check_vs_serial_reference(gi,name):
    from generic import obj
    sr = obj(get_serial_references()[name])
    sg = gi.serial()
    assert(object_eq(sg,sr))
#end def check_vs_serial_reference



def test_files():
    filenames = [
        'rhf.inp',
        'cisd.inp',
        'cas.inp',
        'H.BFD_V5Z_ANO.gms',
        'O.BFD_V5Z.gms',
        ]
    files = get_files()
    assert(set(filenames)==set(files.keys()))
#end def test_files



def test_import():
    from gamess_input import GamessInput,generate_gamess_input
#end def test_import



def test_keyspec_groups():
    from gamess_input import check_keyspec_groups

    divert_nexus_log()

    check_keyspec_groups()

    restore_nexus_log()
#end def test_keyspec_groups



def test_empty_init():
    from gamess_input import GamessInput,generate_gamess_input

    gi = GamessInput()
    assert(len(gi)==0)

    gi = generate_gamess_input()
    assert(isinstance(gi,GamessInput))
    assert(len(gi)==0)
#end def test_empty_init



def test_read():
    from gamess_input import GamessInput
    
    files = get_files()

    input_files = ['rhf.inp','cisd.inp','cas.inp']

    for infile in input_files:
        gi_read = GamessInput(files[infile])
        check_vs_serial_reference(gi_read,infile)
    #end for

#end def test_read



def test_write():
    import os
    from gamess_input import GamessInput

    tpath = testing.setup_unit_test_output_directory('gamess_input','test_write')
    
    files = get_files()

    input_files = ['rhf.inp','cisd.inp','cas.inp']

    for infile in input_files:
        write_file = os.path.join(tpath,infile)

        gi_read = GamessInput(files[infile])

        gi_read.write(write_file)

        gi_write = GamessInput(write_file)

        check_vs_serial_reference(gi_write,infile)
    #end for

#end def test_write



def test_generate():
    import os
    from generic import obj
    from nexus_base import nexus_noncore
    from pseudopotential import Pseudopotentials
    from physical_system import generate_physical_system
    from gamess_input import generate_gamess_input
    
    ppfiles = ['H.BFD_V5Z_ANO.gms','O.BFD_V5Z.gms']

    tpath = testing.setup_unit_test_output_directory(
        test      = 'gamess_input',
        subtest   = 'test_generate',
        divert    = True,
        file_sets = {
            'pseudopotentials':ppfiles
            }
        )

    ppfiles_full = [os.path.join(tpath,'pseudopotentials',f) for f in ppfiles]

    nexus_noncore.pseudopotentials = Pseudopotentials(ppfiles_full)

    input_files = ['rhf.inp','cisd.inp','cas.inp']

    h2o = generate_physical_system(
        elem        = ['O','H','H'], 
        pos         = [[0.000000, 0.000000, 0.000000],
                       [0.000000,-0.757160, 0.586260],
                       [0.000000, 0.757160, 0.586260]],
        units       = 'A',
        net_spin    = 0,  
        O           = 6,  
        H           = 1,  
        # C2v symmetry structure
        folded_elem = ['O','H'],     
        folded_pos  = [[0.000000, 0.000000, 0.000000],
                       [0.000000, 0.757160, 0.586260]],
        )

    inputs = dict()

    inputs['rhf.inp'] = obj(
        system     = h2o,
        pseudos    = ppfiles,
        scftyp     = 'rohf',
        runtyp     = 'energy',
        exetyp     = 'run',
        ispher     = 1,
        maxit      = 200,
        memory     = 150000000,
        dirscf     = True,
        guess      = 'huckel',
        symmetry   = 'Cnv 2',
        )

    inputs['cisd.inp'] = obj(
        system     = h2o,
        pseudos    = ppfiles,
        scftyp     = 'none',
        cityp      = 'guga',
        runtyp     = 'energy',
        exetyp     = 'run',
        ispher     = 1,
        maxit      = 200,
        memory     = 150000000,
        dirscf     = True,
        symmetry   = 'Cnv 2',
        cidrt = obj(
            group  = 'c2v',
            nfzc   = 0,
            ndoc   = 4,
            nalp   = 0,
            nval   = 60,
            nprt   = 2,
            istsym = 1,
            iexcit = 2,
            mxnint = 500000,
            ),
        gugdia = obj(
            prttol = 0.001,
            cvgtol = 1.0e-5,
            itermx = 1000,
            ),
        )

    inputs['cas.inp'] = obj(
        system     = h2o,
        pseudos    = ppfiles,
        scftyp     = 'mcscf',
        runtyp     = 'energy',
        exetyp     = 'run',
        ispher     = 1,
        maxit      = 200,
        memory     = 150000000,
        dirscf     = True,
        symmetry   = 'Cnv 2',
        drt = obj(
            group  = 'c2v',
            nmcc   = 0,
            ndoc   = 4,
            nalp   = 0,
            nval   = 4,
            istsym = 1,
            mxnint = 500000,
            fors   = True,
            ),
        mcscf = obj(
            cistep = 'guga',
            maxit  = 1000,
            fullnr = True,
            acurcy = 1e-5,
            ),
        )

    for infile in input_files:
        gi = generate_gamess_input(**inputs[infile])
        check_vs_serial_reference(gi,infile)
    #end for

    restore_nexus()
#end def test_generate



