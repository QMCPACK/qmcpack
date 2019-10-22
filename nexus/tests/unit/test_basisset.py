
import testing
from testing import value_eq,object_eq
from testing import divert_nexus_log,restore_nexus_log


associated_files = dict()


def get_filenames():
    filenames = [
        'Fe.aug-cc-pwcv5z-dk.0.bas',
        'Fe.aug-cc-pwcv5z-dk.0.gbs',
        'Fe.BFD_VQZ.bas',
        'Fe.BFD_VQZ.gbs',
        'Fe.stuttgart_rsc_1997.0.bas',
        'Fe.stuttgart_rsc_1997.0.gbs',
        'Fe.stuttgart_rsc_1997_ecp.0.bas',
        'Fe.stuttgart_rsc_1997_ecp.0.gbs',
        ]
    return filenames
#end def get_filenames


def get_files():
    return testing.collect_unit_test_file_paths('basisset',associated_files)
#end def get_files



def test_files():
    filenames = get_filenames()
    files = get_files()
    assert(set(files.keys())==set(filenames))
#end def test_files



def test_import():
    import basisset
    from basisset import BasisSets
    from basisset import process_gaussian_text
    from basisset import GaussianBasisSet
#end def test_import



def test_basissets():
    from basisset import BasisSets
    from basisset import BasisFile
    from basisset import gamessBasisFile

    # empty initialization
    BasisSets()
    BasisFile()
    gamessBasisFile()

    filenames = get_filenames()
    files = get_files()
    f = [files[fn] for fn in filenames]

    # standard initialization
    divert_nexus_log()
    bf  = BasisFile(f[1])
    gbf = gamessBasisFile(f[0])
    bs  = BasisSets(f[2:]+[bf,gbf])
    restore_nexus_log()

    assert(bf.element=='Fe')
    assert(bf.filename=='Fe.aug-cc-pwcv5z-dk.0.gbs')

    assert(gbf.element=='Fe')
    assert(gbf.filename=='Fe.aug-cc-pwcv5z-dk.0.bas')
    assert(gbf.text.startswith('IRON'))
    assert(gbf.text.endswith('1         1.3776500              1.0000000'))
    assert(len(gbf.text.strip())==21135)

    for fn in filenames:
        assert(fn in bs)
        assert(isinstance(bs[fn],BasisFile))
        if fn.endswith('.bas'):
            assert(isinstance(bs[fn],gamessBasisFile))
        #end if
    #end for
#end def test_basissets



def test_process_gaussian_text():
    from basisset import process_gaussian_text

    filenames = get_filenames()
    files = get_files()

    bs_ref = {
        'Fe.aug-cc-pwcv5z-dk.0.bas'   : 503, 
        'Fe.aug-cc-pwcv5z-dk.0.gbs'   : 503,
        'Fe.BFD_VQZ.bas'              : 132,
        'Fe.BFD_VQZ.gbs'              : 132,
        'Fe.stuttgart_rsc_1997.0.bas' :  38,
        'Fe.stuttgart_rsc_1997.0.gbs' :  38,

        }

    pp_ref = {
        'Fe.BFD_VQZ.bas'                  : (  9,  132 ),
        'Fe.BFD_VQZ.gbs'                  : ( 13,  132 ),
        'Fe.stuttgart_rsc_1997.0.bas'     : ( 13,   38 ),
        'Fe.stuttgart_rsc_1997.0.gbs'     : ( 17,   38 ),
        'Fe.stuttgart_rsc_1997_ecp.0.bas' : ( 13, None ),
        'Fe.stuttgart_rsc_1997_ecp.0.gbs' : ( 17, None ),
        }

    for fn in filenames:

        if fn.endswith('.bas'):
            format = 'gamess'
        elif fn.endswith('.gbs'):
            format = 'gaussian'
        else:
            format = None
        #end if

        f = open(files[fn],'r')
        text = f.read()
        f.close()

        bs = process_gaussian_text(text,format,pp=False)

        if fn in bs_ref:
            assert(len(bs)==bs_ref[fn])
        #end if

        pp,bs = process_gaussian_text(text,format)

        if fn in pp_ref:
            ppr,bsr = pp_ref[fn]
            assert(len(pp)==ppr)
            if bsr is None:
                assert(bs is None)
            else:
                assert(len(bs)==bsr)
            #end if
        #end if
    #end for
#end def test_process_gaussian_text




def test_gaussianbasisset():
    from basisset import GaussianBasisSet

    filenames = get_filenames()
    files = get_files()

    GaussianBasisSet()

    ref = {
        'Fe.aug-cc-pwcv5z-dk.0.bas'   : 49,
        'Fe.aug-cc-pwcv5z-dk.0.gbs'   : 49,
        'Fe.BFD_VQZ.bas'              : 23,
        'Fe.BFD_VQZ.gbs'              : 23,
        'Fe.stuttgart_rsc_1997.0.bas' : 15,
        'Fe.stuttgart_rsc_1997.0.gbs' : 15,
        }

    gbs = dict()
    for fn in filenames:
        if 'ecp' not in fn:
            if fn.endswith('.bas'):
                format = 'gamess'
            elif fn.endswith('.gbs'):
                format = 'gaussian'
            else:
                format = None
            #end if

            bs = GaussianBasisSet(files[fn],format)

            assert(bs.name=='Fe')
            assert(len(bs.basis)==ref[fn])

            text = bs.write(format=format)
            text = 'header line\n'+text
            bs2 = GaussianBasisSet()
            bs2.read_text(text,format=format)

            assert(object_eq(bs.basis,bs2.basis))
        #end if
    #end for
#end def test_gaussianbasisset

