from ..testing import object_eq
from pathlib import Path


TEST_FILES = [
    Path(__file__+"/../test_basisset_files/Fe.aug-cc-pwcv5z-dk.0.bas").resolve(),
    Path(__file__+"/../test_basisset_files/Fe.aug-cc-pwcv5z-dk.0.gbs").resolve(),
    Path(__file__+"/../test_basisset_files/Fe.BFD_VQZ.bas").resolve(),
    Path(__file__+"/../test_basisset_files/Fe.BFD_VQZ.gbs").resolve(),
    Path(__file__+"/../test_basisset_files/Fe.stuttgart_rsc_1997.0.bas").resolve(),
    Path(__file__+"/../test_basisset_files/Fe.stuttgart_rsc_1997.0.gbs").resolve(),
    Path(__file__+"/../test_basisset_files/Fe.stuttgart_rsc_1997_ecp.0.bas").resolve(),
    Path(__file__+"/../test_basisset_files/Fe.stuttgart_rsc_1997_ecp.0.gbs").resolve(),
]


def test_basissets():
    from ..basisset import BasisSets
    from ..basisset import BasisFile
    from ..basisset import gamessBasisFile

    # empty initialization
    BasisSets()
    BasisFile()
    gamessBasisFile()

    basis_file  = BasisFile(TEST_FILES[1])
    gamess_basis_file = gamessBasisFile(TEST_FILES[0])
    basis_sets  = BasisSets(TEST_FILES[2:]+[basis_file,gamess_basis_file])

    assert(basis_file.element=='Fe')
    assert(basis_file.filename=='Fe.aug-cc-pwcv5z-dk.0.gbs')

    assert(gamess_basis_file.element=='Fe')
    assert(gamess_basis_file.filename=='Fe.aug-cc-pwcv5z-dk.0.bas')
    assert(gamess_basis_file.text.startswith('IRON'))
    assert(gamess_basis_file.text.endswith('1         1.3776500              1.0000000'))
    assert(len(gamess_basis_file.text.strip())==21135)

    for fn in TEST_FILES:
        assert(fn.name in basis_sets)
        assert(isinstance(basis_sets[fn.name],BasisFile))
        if fn.suffix == '.bas':
            assert(isinstance(basis_sets[fn.name],gamessBasisFile))
        #end if
    #end for
#end def test_basissets



def test_process_gaussian_text():
    from ..basisset import process_gaussian_text

    basis_refs = {
        'Fe.aug-cc-pwcv5z-dk.0.bas'   : 503, 
        'Fe.aug-cc-pwcv5z-dk.0.gbs'   : 503,
        'Fe.BFD_VQZ.bas'              : 132,
        'Fe.BFD_VQZ.gbs'              : 132,
        'Fe.stuttgart_rsc_1997.0.bas' :  38,
        'Fe.stuttgart_rsc_1997.0.gbs' :  38,

        }

    pseudo_refs = {
        'Fe.BFD_VQZ.bas'                  : (  9,  132 ),
        'Fe.BFD_VQZ.gbs'                  : ( 13,  132 ),
        'Fe.stuttgart_rsc_1997.0.bas'     : ( 13,   38 ),
        'Fe.stuttgart_rsc_1997.0.gbs'     : ( 17,   38 ),
        'Fe.stuttgart_rsc_1997_ecp.0.bas' : ( 13, None ),
        'Fe.stuttgart_rsc_1997_ecp.0.gbs' : ( 17, None ),
        }

    for file in TEST_FILES:

        if file.suffix == '.bas':
            format = 'gamess'
        elif file.suffix == '.gbs':
            format = 'gaussian'
        else:
            format = None
        #end if

        file_text = file.read_text()

        basis_set = process_gaussian_text(file_text,format,pp=False)

        if file.name in basis_refs:
            assert(len(basis_set)==basis_refs[file.name])
        #end if

        pseudo, basis_set = process_gaussian_text(file_text,format)

        if file.name in pseudo_refs:
            pseudo_ref, basis_ref = pseudo_refs[file.name]
            assert(len(pseudo)==pseudo_ref)
            if basis_ref is not None:
                assert(len(basis_set)==basis_ref)
            #end if
        #end if
    #end for
#end def test_process_gaussian_text




def test_gaussianbasisset():
    from ..basisset import GaussianBasisSet

    GaussianBasisSet()

    ref = {
        'Fe.aug-cc-pwcv5z-dk.0.bas'   : 49,
        'Fe.aug-cc-pwcv5z-dk.0.gbs'   : 49,
        'Fe.BFD_VQZ.bas'              : 23,
        'Fe.BFD_VQZ.gbs'              : 23,
        'Fe.stuttgart_rsc_1997.0.bas' : 15,
        'Fe.stuttgart_rsc_1997.0.gbs' : 15,
        }

    for file in TEST_FILES:
        if 'ecp' not in file.name:
            if file.suffix == '.bas':
                format = 'gamess'
            elif file.suffix == '.gbs':
                format = 'gaussian'
            else:
                format = None
            #end if

            bs = GaussianBasisSet(file,format)

            assert(bs.name=='Fe')
            assert(len(bs.basis)==ref[file.name])

            text = bs.write(format=format)
            text = 'header line\n'+text
            bs2 = GaussianBasisSet()
            bs2.read_text(text,format=format)

            assert(object_eq(bs.basis,bs2.basis))
        #end if
    #end for
#end def test_gaussianbasisset

