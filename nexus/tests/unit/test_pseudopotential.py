
import testing
from testing import divert_nexus_log,restore_nexus_log


associated_files = dict()


def get_filenames():
    filenames = [
        'C.BFD.gms',
        'C.BFD.upf',
        'C.BFD.xml',
        ]
    return filenames
#end def get_filenames


def get_files():
    return testing.collect_unit_test_file_paths('pseudopotential',associated_files)
#end def get_files



def test_files():
    filenames = get_filenames()
    files = get_files()
    assert(set(files.keys())==set(filenames))
#end def test_files


def test_import():
    import pseudopotential
    from pseudopotential import Pseudopotentials
    from pseudopotential import PseudoFile
    from pseudopotential import gamessPPFile
    from pseudopotential import PPset
    from pseudopotential import Pseudopotential
    from pseudopotential import SemilocalPP
    from pseudopotential import GaussianPP
    from pseudopotential import QmcpackPP
    from pseudopotential import CasinoPP
#end def test_import



def test_pp_elem_label():
    from pseudopotential import pp_elem_label

    ppfiles = [
        ('C'  , 'C'  , 'C.BFD.xml'   ),
        ('C2' , 'C'  , 'C2.ONCV.upf' ),
        ('C'  , 'C'  , 'C_ONCV.upf'  ),
        ('C'  , 'C'  , 'C-ONCV.upf'  ),
        ('Co' , 'Co' , 'Co.BFD.xml'  ),
        ('Co2', 'Co' , 'Co2.ONCV.upf'),
        ('Co' , 'Co' , 'Co_ONCV.upf' ),
        ('Co' , 'Co' , 'Co-ONCV.upf' ),
        ]

    for reflabel,refsymbol,ppfile in ppfiles:
        label,symbol,is_elem = pp_elem_label(ppfile)
        assert(is_elem)
        assert(label==reflabel)
        assert(symbol==refsymbol)
    #end for

#end def test_pp_elem_label



def test_pseudopotentials():
    from pseudopotential import Pseudopotentials
    from pseudopotential import PseudoFile
    from pseudopotential import gamessPPFile

    filenames = get_filenames()
    files = get_files()

    filepaths = [files[fn] for fn in filenames]

    # empty initialization
    Pseudopotentials()
    PseudoFile()
    gamessPPFile()

    # standard initialization
    divert_nexus_log()
    pps = Pseudopotentials(filepaths)
    restore_nexus_log()

    for fn in filenames:
        assert(fn in pps)
        pp = pps[fn]
        assert(isinstance(pp,PseudoFile))
        if fn.endswith('.gms'):
            assert(isinstance(pp,gamessPPFile))
        #end if
        assert(pp.element=='C')
        assert(pp.element_label=='C')
        assert(pp.filename==fn)
    #end for

    basis_ref = '''s 9 1.00
1 0.051344     0.013991
2 0.102619     0.169852
3 0.205100     0.397529
4 0.409924     0.380369
5 0.819297     0.180113
6 1.637494     -0.033512
7 3.272791     -0.121499
8 6.541187     0.015176
9 13.073594     -0.000705
s 1 1.00
1 0.098302     1.000000
s 1 1.00
1 0.232034     1.000000
s 1 1.00
1 0.744448     1.000000
s 1 1.00
1 1.009914     1.000000
p 9 1.00
1 0.029281     0.001787
2 0.058547     0.050426
3 0.117063     0.191634
4 0.234064     0.302667
5 0.468003     0.289868
6 0.935757     0.210979
7 1.871016     0.112024
8 3.741035     0.054425
9 7.480076     0.021931
p 1 1.00
1 0.084047     1.000000
p 1 1.00
1 0.216618     1.000000
p 1 1.00
1 0.576869     1.000000
p 1 1.00
1 1.006252     1.000000
d 1 1.00
1 0.206619     1.000000
d 1 1.00
1 0.606933     1.000000
d 1 1.00
1 1.001526     1.000000
d 1 1.00
1 1.504882     1.000000
f 1 1.00
1 0.400573     1.000000
f 1 1.00
1 1.099564     1.000000
f 1 1.00
1 1.501091     1.000000
g 1 1.00
1 0.797648     1.000000
g 1 1.00
1 1.401343     1.000000
h 1 1.00
1 1.001703     1.000000'''

    pp_ref = '''C-QMC GEN 2 1
3
4.00000000 1 8.35973821
33.43895285 3 4.48361888
-19.17537323 2 3.93831258
1
22.55164191 2 5.02991637'''

    pp = pps['C.BFD.gms']
    assert(pp.basis_text==basis_ref)
    assert(pp.pp_text==pp_ref)

#end def test_pseudopotentials
