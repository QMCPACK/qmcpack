import pytest
from . import NexusTestOrder
pytestmark = pytest.mark.order(NexusTestOrder.VASP_ANALYZER)


from ..generic import NexusError
from ..pseudoset import PseudoSet

from nexus.nexus_base import nexus_core
from . import isolate_nexus_core, TEST_DIR
from .. import testing
from ..testing import object_eq,dict_serialize


serial_references = dict()

c_potcar_text = '''
This is not a real POTCAR file.

End of Dataset
'''

def generate_serial_references():
    import numpy as np
    from ..developer import obj

    ref = {
        'incar/encut' : 450.0,
        'incar/ibrion' : 2,
        'incar/icharg' : 2,
        'incar/isif' : 2,
        'incar/istart' : 0,
        'incar/nelect' : 64,
        'incar/nsw' : 5,
        'kpoints/centering' : 'monkhorst',
        'kpoints/kgrid' : (2, 2, 2),
        'kpoints/kshift' : (0, 0, 0),
        'kpoints/mode' : 'auto',
        'poscar/axes' : np.array([
            [3.57, 3.57, 0.0],
            [0.0, 3.57, 3.57],
            [3.57, 0.0, 3.57]]),
        'poscar/coord' : 'cartesian',
        'poscar/description' : None,
        'poscar/dynamic' : None,
        'poscar/elem' : ['C'],
        'poscar/elem_count' : [16],
        'poscar/pos' : np.array([
            [0.0, 0.0, 0.0],
            [0.8925, 0.8925, 0.8925],
            [1.785, 1.785, 0.0],
            [2.6775, 2.6775, 0.8925],
            [0.0, 1.785, 1.785],
            [0.8925, 2.6775, 2.6775],
            [1.785, 3.57, 1.785],
            [2.6775, 4.4625, 2.6775],
            [1.785, 0.0, 1.785],
            [2.6775, 0.8925, 2.6775],
            [3.57, 1.785, 1.785],
            [4.4625, 2.6775, 2.6775],
            [1.785, 1.785, 3.57],
            [2.6775, 2.6775, 4.4625],
            [3.57, 3.57, 3.57],
            [4.4625, 4.4625, 4.4625]]),
        'poscar/scale' : 1.0,
        'poscar/vel' : None,
        'poscar/vel_coord' : None,
        'potcar/files' : ['C.POTCAR'],
        'potcar/pseudos' : obj(),
        }

    ref_read = ref.copy()

    ref_read['incar/nelect']       = 64.0
    ref_read['kpoints/centering']  = 'monkhorst-pack'
    ref_read['kpoints/kgrid']      = np.array((2, 2, 2),dtype=int)
    ref_read['kpoints/kshift']     = np.array((0, 0, 0),dtype=float)
    ref_read['poscar/description'] = 'System cell and coordinates'
    ref_read['poscar/elem']        = np.array(['C'],dtype=str)
    ref_read['poscar/elem_count']  = np.array([16],dtype=int)
    ref_read['potcar/files']       = None
    del ref_read['potcar/pseudos']
    ref_read['potcar/pseudos/0']   = c_potcar_text

    serial_references['read']     = ref_read

    serial_references['write']    = ref_read.copy()

    serial_references['generate'] = ref.copy()

#end def generate_serial_references


def get_serial_references():
    if len(serial_references)==0:
        generate_serial_references()
    #end if
    return serial_references
#end def get_serial_references


def check_vs_serial_reference(gi,name):
    from ..developer import obj
    sr = obj(get_serial_references()[name])
    sg = dict_serialize(gi,dict_type=obj)
    same = object_eq(sg,sr)
    if not same:
        print('\n'+name+' differs')
        testing.print_diff(sr,sg)
    #end if
    assert(same)
#end def check_vs_serial_reference


TEST_FILES = {
    "d16bulk.POSCAR":  TEST_DIR / "test_vasp_input_files/d16bulk.POSCAR",
    "diamond_INCAR":   TEST_DIR / "test_vasp_input_files/diamond_INCAR",
    "diamond_KPOINTS": TEST_DIR / "test_vasp_input_files/diamond_KPOINTS",
    "diamond_POSCAR":  TEST_DIR / "test_vasp_input_files/diamond_POSCAR",
    "diamond_POTCAR":  TEST_DIR / "test_vasp_input_files/diamond_POTCAR",
    }

for file in TEST_FILES.values():
    assert(file.exists()), f"Test file not found! {file}"


def test_keyword_consistency():
    from ..vasp_input import Incar,Stopcar,VKeywordFile

    for cls in Incar,Stopcar:
        assert(isinstance(cls.scalar_keywords,frozenset))
        assert(isinstance(cls.array_keywords,frozenset))
        assert(isinstance(cls.keywords,frozenset))
        scalar_keywords = set()
        array_keywords = set()
        typed_keywords = set()
        for keyword_type in cls.kw_scalars+cls.kw_arrays:
            keywords = getattr(cls,keyword_type)
            assert(typed_keywords.isdisjoint(keywords))
            typed_keywords |= keywords
            if keyword_type in cls.kw_scalars:
                scalar_keywords |= keywords
            else:
                array_keywords |= keywords
            #end if
        #end for
        assert(scalar_keywords==cls.scalar_keywords)
        assert(array_keywords==cls.array_keywords)
        assert(typed_keywords==cls.keywords)
        assert(set(cls.dtype.keys())==cls.keywords)
        assert(set(cls.read_value.keys())==cls.keywords)
        assert(set(cls.write_value.keys())==cls.keywords)
        assert(set(cls.assign_value.keys())==cls.keywords)
    #end for

    class AdditionalKeywordFile(VKeywordFile):
        bools = frozenset({'flag'})
    #end class AdditionalKeywordFile

    AdditionalKeywordFile.class_init()
    assert(AdditionalKeywordFile.keywords==frozenset({'flag'}))
    assert(AdditionalKeywordFile.dtype.flag=='bools')
#end def test_keyword_consistency


def test_current_keyword_types():
    from ..vasp_input import Incar

    expected = {
        'ints': {
            'efermi_nedos', 'elph_nbands', 'iopt', 'libmbd_n_omega_grid',
            'ml_output_mode',
            },
        'reals': {
            'apaco', 'ddr', 'elph_kspacing', 'nupdown', 'spring',
            },
        'bools': {
            'elph_run', 'lglobal', 'ml_lib',
            },
        'strings': {
            'elph_driver', 'ml_mode', 'pyamff_model', 'xc',
            },
        'int_arrays': {
            'elph_nbands_sum', 'libmbd_k_grid',
            },
        'real_arrays': {
            'elph_selfen_temps', 'pomass', 'quad_efg', 'smearings',
            },
        'bool_arrays': {
            'fmp_active',
            },
        }

    for keyword_type,keywords in expected.items():
        for keyword in keywords:
            assert(Incar.dtype[keyword]==keyword_type)
        #end for
    #end for

    expected_mixed = {
        'bext':   ('reals','real_arrays'),
        'efermi': ('strings','reals'),
        'libxc1': ('strings','ints'),
        'libxc2': ('strings','ints'),
        'lreal':  ('strings','bools'),
        }
    assert(dict(Incar.mixed_types)==expected_mixed)
    for keyword,types in expected_mixed.items():
        assert(Incar.dtype[keyword]==types)
    #end for
#end def test_current_keyword_types


def test_current_keywords_roundtrip():
    import numpy as np
    from ..vasp_input import Incar

    incar = Incar()
    incar.assign(
        elph_driver       = 'el',
        elph_nbands       = 8,
        elph_run          = True,
        elph_selfen_temps = [0.0,300.0],
        fmp_active        = [True,False],
        iopt              = 7,
        libmbd_k_grid     = [2,2,1],
        lglobal           = True,
        spring            = -5.0,
        )

    reread = Incar()
    reread.read_text(incar.write_text())

    assert(reread.elph_driver=='el')
    assert(reread.elph_nbands==8)
    assert(reread.elph_run)
    assert(np.array_equal(reread.elph_selfen_temps,[0.0,300.0]))
    assert(np.array_equal(reread.fmp_active,[True,False]))
    assert(reread.iopt==7)
    assert(np.array_equal(reread.libmbd_k_grid,[2,2,1]))
    assert(reread.lglobal)
    assert(reread.spring==-5.0)
#end def test_current_keywords_roundtrip


def test_mixed_keyword_types():
    import numpy as np
    from ..vasp_input import (
        Incar,assign_mixed,generate_vasp_input,mixed_type_matches,
        )

    generated = generate_vasp_input(
        bext   = 0.02,
        efermi = 5.4,
        libxc1 = 101,
        libxc2 = 'GGA_C_PBE',
        lreal  = False,
        )
    incar = generated.incar

    assert(incar.bext==0.02)
    assert(incar.efermi==5.4)
    assert(incar.libxc1==101)
    assert(incar.libxc2=='GGA_C_PBE')
    assert(incar.lreal is False)

    reread = Incar()
    reread.read_text(incar.write_text())
    assert(reread.bext==0.02)
    assert(reread.efermi==5.4)
    assert(reread.libxc1==101)
    assert(reread.libxc2=='GGA_C_PBE')
    assert(reread.lreal is False)

    reread.read_text(
        'BEXT = 0 0 0.1\n'
        'EFERMI = MIDGAP\n'
        'LIBXC1 = GGA_X_PBE\n'
        'LREAL = Auto\n'
        )
    assert(np.array_equal(reread.bext,[0.0,0.0,0.1]))
    assert(reread.efermi=='MIDGAP')
    assert(reread.libxc1=='GGA_X_PBE')
    assert(reread.lreal=='Auto')

    invalid = (
        ('bext','0.1'),
        ('efermi',True),
        ('libxc1',1.25),
        ('lreal',2),
        )
    for keyword,value in invalid:
        with pytest.raises(
            NexusError,
            match='assign failed for keyword {0}'.format(keyword),
            ):
            Incar().assign(**{keyword:value})
        #end with
    #end for

    with pytest.raises(NexusError,match='read failed for keyword efermi'):
        Incar().read_text('EFERMI = 1.0 2.0')
    #end with

    with pytest.raises(
        ValueError,
        match=r'unknown value for keyword value_type.*must be one of',
        ):
        mixed_type_matches(1.0,'unknown')
    #end with
    with pytest.raises(ValueError,match='ints for value 1.5'):
        assign_mixed(1.5,types=('ints',))
    #end with
#end def test_mixed_keyword_types


def test_block_constructs():
    from types import MappingProxyType
    from ..developer import obj
    from ..vasp_input import Incar,generate_vasp_input

    assert(isinstance(Incar.block_constructs,MappingProxyType))
    assert(isinstance(
        Incar.block_constructs['kernel_truncation'],MappingProxyType
        ))
    assert(Incar.block_constructs['kernel_truncation']['lcoarsen']=='bool')
    assert(Incar.block_constructs['plugins']['ml_outblock']=='int')
    assert(Incar.block_constructs['image_*']=='keywords')

    generated = generate_vasp_input(
        kernel_truncation = obj(
            factor          = 0.5,
            idimensionality = 2,
            isurface        = 3,
            lcoarsen        = False,
            ltruncate       = True,
            ),
        plugins = {
            'force_and_stress': True,
            'ml_mode':          'run',
            'ml_outblock':      4,
            'neighbor_cutoff':  6.5,
            },
        )
    incar = generated.incar
    assert(isinstance(incar.kernel_truncation,obj))
    assert(isinstance(incar.plugins,obj))

    text = incar.write_text()
    assert(
        'KERNEL_TRUNCATION = {\n'
        '  FACTOR          = 0.5\n'
        '  IDIMENSIONALITY = 2\n'
        '  ISURFACE        = 3\n'
        '  LCOARSEN        = .FALSE.\n'
        '  LTRUNCATE       = .TRUE.\n'
        '}\n' in text
        )
    assert(
        'PLUGINS = {\n'
        '  FORCE_AND_STRESS = .TRUE.\n'
        '  ML_MODE          = run\n'
        '  ML_OUTBLOCK      = 4\n'
        '  NEIGHBOR_CUTOFF  = 6.5\n'
        '}\n' in text
        )

    reread = Incar()
    reread.read_text(text)
    assert(isinstance(reread.kernel_truncation,obj))
    assert(reread.kernel_truncation==incar.kernel_truncation)
    assert(reread.plugins==incar.plugins)

    with pytest.raises(
        TypeError,
        match='block construct value should be a mapping, but is list',
        ):
        Incar().assign_block_construct('plugins',[])
    #end with
#end def test_block_constructs


def test_block_construct_parsing():
    import numpy as np
    from ..developer import obj
    from ..vasp_input import Incar

    text = (
        'ENCUT = 400\n'
        'KERNEL_TRUNCATION {\n'
        '  LTRUNCATE = T; IDIMENSIONALITY = 2\n'
        '  ISURFACE = 3 # surface normal\n'
        '}\n'
        'PLUGINS/FORCE_AND_STRESS = T\n'
        'PLUGINS/ML_MODE = "run"\n'
        'IMAGE_2 = {\n'
        '  ENCUT = 450\n'
        '  LREAL = Auto\n'
        '  MAGMOM = 2*1D-1\n'
        '}\n'
        )
    incar = Incar()
    incar.read_text(text)

    assert(incar.encut==400.0)
    assert(incar.kernel_truncation==obj(
        ltruncate=True,idimensionality=2,isurface=3
        ))
    assert(incar.plugins==obj(force_and_stress=True,ml_mode='run'))
    assert(isinstance(incar.image_2,obj))
    assert(incar.image_2.encut==450.0)
    assert(incar.image_2.lreal=='Auto')
    assert(np.array_equal(incar.image_2.magmom,[0.1,0.1]))

    reread = Incar()
    reread.read_text(incar.write_text())
    assert(reread.kernel_truncation==incar.kernel_truncation)
    assert(reread.plugins==incar.plugins)
    assert(reread.image_2.encut==450.0)
    assert(np.array_equal(reread.image_2.magmom,[0.1,0.1]))

    ordered = Incar()
    ordered.read_text(
        'PLUGINS/FORCE_AND_STRESS = T\n'
        'PLUGINS/ML_MODE = none\n'
        'PLUGINS = { ML_MODE = run }\n'
        )
    assert(ordered.plugins.force_and_stress)
    assert(ordered.plugins.ml_mode=='run')

    with pytest.raises(ValueError,match='not a field of block construct'):
        Incar().assign(kernel_truncation={'unknown':1})
    #end with
    with pytest.raises(ValueError,match='invalid literal for int'):
        Incar().read_text('PLUGINS = { ML_OUTBLOCK = 1.5 }')
    #end with
    with pytest.raises(
        NexusError,
        match='block construct plugins is not closed',
        ):
        Incar().read_text('PLUGINS = { ML_MODE = run')
    #end with
#end def test_block_construct_parsing


def test_integer_assignment_validation():
    import numpy as np
    from ..vasp_input import Incar,assign_int_array,generate_vasp_input

    incar = Incar()
    incar.assign(
        ibrion      = np.int64(2),
        images      = 3.0+5e-9,
        random_seed = [1,2.0-5e-9,np.int32(3)],
        )
    assert(incar.ibrion==2)
    assert(incar.images==3)
    assert(np.array_equal(incar.random_seed,[1,2,3]))

    generated = generate_vasp_input(ibrion=2.0+5e-9)
    assert(generated.incar.ibrion==2)

    invalid_scalars = (
        True,
        '3',
        2.0+2e-8,
        2.1,
        np.nan,
        np.inf,
        float(2**54),
        )
    for value in invalid_scalars:
        with pytest.raises(NexusError,match='assign failed for keyword ibrion'):
            Incar().assign(ibrion=value)
        #end with
    #end for

    invalid_arrays = (
        [1,2.1,3],
        [1,True,3],
        [1,'2',3],
        [1,2**100,3],
        [[1,2],[3,4]],
        )
    for value in invalid_arrays:
        with pytest.raises(
            NexusError,
            match='assign failed for keyword random_seed',
            ):
            Incar().assign(random_seed=value)
        #end with
    #end for

    with pytest.raises(NexusError,match='element at index 1'):
        Incar().assign(random_seed=[1,2.1,3])
    #end with
    with pytest.raises(OverflowError):
        assign_int_array([1,2**100,3])
    #end with
    with pytest.raises(NexusError,match='assign failed for keyword ibrion'):
        generate_vasp_input(ibrion=2.1)
    #end with
#end def test_integer_assignment_validation


def test_keyword_syntax_parsing():
    import numpy as np
    from ..vasp_input import Incar

    text = (
        'ENCUT = 4D2; ISMEAR = -1 ! trailing comment\n'
        'MAGMOM = 2*1.5 \\\n'
        '         2*-1.5 # another comment\n'
        'LCHARG = .false.\n'
        'IBRION = 1; IBRION = 2\n'
        )
    incar = Incar()
    incar.read_text(text)

    assert(incar.encut==400.0)
    assert(incar.ismear==-1)
    assert(np.array_equal(incar.magmom,[1.5,1.5,-1.5,-1.5]))
    assert(not incar.lcharg)
    assert(incar.ibrion==2)
#end def test_keyword_syntax_parsing


def test_quoted_string_preservation():
    from ..vasp_input import Incar

    text = (
        'SYSTEM = "a; b # c ! d = e"; '
        'WANNIER90_WIN = "Begin; x=y # z!"\n'
        'ENCUT = 400\n'
        )
    incar = Incar()
    incar.read_text(text)

    assert(incar.system=='a; b # c ! d = e')
    assert(incar.wannier90_win=='Begin; x=y # z!')
    assert(incar.encut==400.0)

    generated = Incar()
    generated.assign(
        system        = ' a; b # c ! d ',
        wannier90_win = 'line ending in \\',
        )
    generated_text = generated.write_text()

    assert('SYSTEM        = " a; b # c ! d "' in generated_text)
    assert('WANNIER90_WIN = "line ending in \\"' in generated_text)

    reread = Incar()
    reread.read_text(generated_text)
    assert(reread.system==' a; b # c ! d ')
    assert(reread.wannier90_win=='line ending in \\')
#end def test_quoted_string_preservation


def test_array_compression_roundtrip():
    import numpy as np
    from ..vasp_input import Incar

    incar = Incar()
    incar.assign(
        lattice_constraints = [True]*4,
        magmom               = [1.5]*4,
        random_seed          = [2]*4,
        )
    text = incar.write_text()

    assert('4*T' in text)
    assert('4*1.5' in text)
    assert('4*2' in text)

    reread = Incar()
    reread.read_text(text)
    assert(np.array_equal(reread.lattice_constraints,[True]*4))
    assert(np.array_equal(reread.magmom,[1.5]*4))
    assert(np.array_equal(reread.random_seed,[2]*4))

    values = np.array([1.0,1.0000004,1.0000002,2.0,2.0,2.0,2.0])
    incar.assign(magmom=values)
    text = incar.write_text()

    assert('MAGMOM              = 1.0 1.0000004 1.0000002 4*2.0' in text)

    reread.read_text(text)
    assert(np.array_equal(reread.magmom,values))
#end def test_array_compression_roundtrip


def test_real_array_fortran_exponents():
    import numpy as np
    from ..vasp_input import Incar

    incar = Incar()
    incar.read_text('MAGMOM = 2*1D-1 -2d+0')

    assert(np.array_equal(incar.magmom,[0.1,0.1,-2.0]))
#end def test_real_array_fortran_exponents


def test_invalid_keyword_input():
    from ..vasp_input import Incar

    with pytest.raises(NexusError,match='NOT_A_TAG is not a keyword'):
        Incar().assign(not_a_tag=1)
    #end with
    with pytest.raises(NexusError,match='read failed for keyword lcharg'):
        Incar().read_text('LCHARG = maybe')
    #end with
    with pytest.raises(NexusError,match='assign failed for keyword magmom'):
        Incar().assign(magmom=1.0)
    #end with
    with pytest.raises(NexusError,match=r'quotation marks.*not paired'):
        Incar().read_text('SYSTEM = "unfinished')
    #end with
    with pytest.raises(NexusError,match='incomplete line continuation'):
        Incar().read_text('ENCUT = 400 \\')
    #end with
    with pytest.raises(NexusError,match='repeat count must be non-negative'):
        Incar().read_text('MAGMOM = -2*1.0')
    #end with
#end def test_invalid_keyword_input


def test_keyword_catalog_integrity():
    from ..vasp_input import Incar

    keyword_sets = (
        Incar.ints,
        Incar.reals,
        Incar.bools,
        Incar.strings,
        Incar.int_arrays,
        Incar.real_arrays,
        Incar.bool_arrays,
        )
    assert(len(Incar.keywords)==687)
    assert(sum(map(len,keyword_sets))==len(Incar.keywords))
    assert(Incar.deprecated<=Incar.keywords)
    assert(Incar.unsupported.isdisjoint(Incar.keywords))
    assert('elph_transport_nedos' in Incar.ints)
    assert('csvr_period' in Incar.reals)
    assert('elph_rotateprojectors' in Incar.bools)
#end def test_keyword_catalog_integrity


def test_deprecated_keyword_compatibility():
    from ..vasp_input import Incar

    incar = Incar()
    incar.assign(
        ichain     = 0,
        jacobian   = 1.0,
        lclimb     = True,
        ldneb      = False,
        timestep   = 0.1,
        )
    reread = Incar()
    reread.read_text(incar.write_text())

    assert(reread.ichain==0)
    assert(reread.jacobian==1.0)
    assert(reread.lclimb)
    assert(not reread.ldneb)
    assert(reread.timestep==0.1)
#end def test_deprecated_keyword_compatibility


def test_keyword_file_roundtrip():
    import numpy as np
    from ..vasp_input import Incar

    text = 'SYSTEM = "first line\nsecond line"\nENCUT = 400\nLATTICE_CONSTRAINTS = t f .TRUE.\n'
    incar = Incar()
    incar.read_text(text)

    assert(incar.system=='first line\nsecond line')
    assert(incar.encut==400.0)
    assert(np.array_equal(incar.lattice_constraints,(True,False,True)))

    reread = Incar()
    reread.read_text(incar.write_text())
    assert(reread.system==incar.system)
    assert(reread.encut==incar.encut)
#end def test_keyword_file_roundtrip


def test_kpoints_tetrahedra_roundtrip():
    import numpy as np
    from ..developer import obj
    from ..vasp_input import Kpoints

    kpoints = Kpoints()
    kpoints.mode = 'explicit'
    kpoints.coord = 'reciprocal'
    kpoints.kpoints = np.array([[0.0,0.0,0.0]])
    kpoints.kweights = np.array([1.0])
    kpoints.tetrahedra = obj()
    kpoints.tetrahedra[0] = obj(
        volume     = 0.5,
        degeneracy = 1,
        corners    = np.array([1,1,1,1]),
        )

    reread = Kpoints()
    reread.read_text(kpoints.write_text())

    assert(reread.mode=='explicit')
    assert(len(reread.tetrahedra)==1)
    assert(reread.tetrahedra[0].volume==0.5)
    assert(np.array_equal(reread.tetrahedra[0].corners,(1,1,1,1)))
#end def test_kpoints_tetrahedra_roundtrip


def test_kpoints_labels_and_related_files():
    import numpy as np
    from ..vasp_input import (
        Kpoints,KpointsElph,KpointsOpt,KpointsWan,Qpoints
        )

    text = (
        'band path\n'
        '20\n'
        'line mode\n'
        'reciprocal\n'
        '0 0 0 Gamma\n'
        '0.5 0 0 X\n'
        )
    for cls in (Kpoints,KpointsElph,KpointsOpt,KpointsWan,Qpoints):
        kpoints = cls()
        kpoints.read_text(text)
        assert(np.array_equal(kpoints.kendpoints,[[0,0,0],[0.5,0,0]]))
        assert(np.array_equal(kpoints.labels,['Gamma','X']))

        reread = cls()
        reread.read_text(kpoints.write_text())
        assert(np.array_equal(reread.labels,['Gamma','X']))
    #end for

    automatic = Kpoints()
    automatic.read_text('automatic\n0\nAuto\n10.5\n')
    assert(automatic.kgrid==10.5)
    reread = Kpoints()
    reread.read_text(automatic.write_text())
    assert(reread.kgrid==10.5)

    explicit = Kpoints()
    explicit.read_text(
        'labeled points\n'
        '2\n'
        'reciprocal\n'
        '0 0 0 1 Gamma\n'
        '0.5 0 0 0 X point\n'
        )
    assert(np.array_equal(explicit.labels,['Gamma','X point']))
    reread = Kpoints()
    reread.read_text(explicit.write_text())
    assert(np.array_equal(reread.labels,explicit.labels))

    blank_description = Kpoints()
    blank_description.read_text('\n0\nGamma\n12 12 12\n')
    assert(blank_description.mode=='auto')
    assert(np.array_equal(blank_description.kgrid,[12,12,12]))
#end def test_kpoints_labels_and_related_files


def test_poscar_coordinates_and_structure_conversion():
    import numpy as np
    from ..vasp_input import Poscar,VaspInput

    poscar = Poscar()
    poscar.scale = 2.0
    poscar.axes = np.eye(3)
    poscar.elem = np.array(['H'])
    poscar.elem_count = np.array([1])
    poscar.coord = 'direct'
    poscar.pos = np.array([[0.25,0.25,0.25]])

    vasp_input = VaspInput()
    vasp_input.poscar = poscar
    poscar.change_specifier('cartesian',vasp_input)
    poscar.change_specifier('direct',vasp_input)
    assert(np.allclose(poscar.pos,[[0.25,0.25,0.25]]))
    assert(poscar.check_complete()=='')

    structure = vasp_input.return_system(structure_only=True)
    assert(np.allclose(structure.axes,2.0*np.eye(3)))
    assert(np.allclose(structure.pos,[[0.5,0.5,0.5]]))

    selective_text = '''selective dynamics
1.0
1 0 0
0 1 0
0 0 1
H
1
Selective dynamics
Direct
0 0 0 t f t
'''
    selective = Poscar()
    selective.read_text(selective_text)
    assert(np.array_equal(selective.dynamic,[[True,False,True]]))
#end def test_poscar_coordinates_and_structure_conversion


def test_poscar_extended_format():
    import numpy as np
    from ..vasp_input import Poscar,VaspInput

    text = (
        'extended structure # retained as a comment\n'
        '1 2 3\n'
        '1 0 0\n'
        '0 1 0\n'
        '0 0 1\n'
        'H\n'
        '1\n'
        'direct\n'
        '0.25 0.25 0.25 center\n'
        'Lattice velocities and vectors\n'
        '1\n'
        '0.1 0 0\n'
        '0 0.2 0\n'
        '0 0 0.3\n'
        '1 0 0\n'
        '0 2 0\n'
        '0 0 3\n'
        '\n'
        '1.0 2.0 3.0\n'
        '\n'
        'predictor-corrector state\n'
        '0.5\n'
        )
    poscar = Poscar()
    poscar.read_text(text)

    assert(np.array_equal(poscar.scale,[1.0,2.0,3.0]))
    assert(np.array_equal(poscar.labels,['center']))
    assert(poscar.lattice_vel_init==1)
    assert(np.array_equal(poscar.lattice_vel,np.diag([0.1,0.2,0.3])))
    assert(poscar.vel_coord=='cartesian')
    assert(np.array_equal(poscar.vel,[[1.0,2.0,3.0]]))
    assert('predictor-corrector state' in poscar.md_extra)

    vasp_input = VaspInput()
    vasp_input.poscar = poscar
    structure = vasp_input.return_system(structure_only=True)
    assert(np.array_equal(structure.axes,np.diag([1.0,2.0,3.0])))
    assert(np.allclose(structure.pos,[[0.25,0.5,0.75]]))
    assert(np.array_equal(structure.vel,[[1.0,2.0,3.0]]))

    reread = Poscar()
    reread.read_text(poscar.write_text())
    assert(np.array_equal(reread.scale,poscar.scale))
    assert(np.array_equal(reread.labels,poscar.labels))
    assert(reread.vel_coord=='cartesian')
    assert(reread.md_extra==poscar.md_extra)

    poscar.elem = np.array(['H1'])
    assert('H1 ' in poscar.write_text())
#end def test_poscar_extended_format


def test_auxiliary_formatted_files():
    import numpy as np
    from ..vasp_input import Iconst,Ircar,Irccar,Penaltypot

    assert(Ircar is Irccar)

    iconst = Iconst()
    iconst.read_text(
        'R 1 2 0\n'
        'W 1 2 1.1 9 14 5\n'
        'S 1.0 -1.0 7\n'
        )
    assert(iconst.coordinates[0].flag=='R')
    assert(iconst.coordinates[1]['items']==(1,2,1.1,9,14))
    assert(iconst.coordinates[2].status==7)
    reread_iconst = Iconst()
    reread_iconst.read_text(iconst.write_text())
    assert(reread_iconst.coordinates==iconst.coordinates)

    penaltypot = Penaltypot()
    penaltypot.read_text('1.0 2.0 0.5\n2.0 3.0 -0.5\n')
    assert(np.array_equal(
        penaltypot.hills,[[1.0,2.0,0.5],[2.0,3.0,-0.5]]
        ))
    reread_penalty = Penaltypot()
    reread_penalty.read_text(penaltypot.write_text())
    assert(np.array_equal(reread_penalty.hills,penaltypot.hills))

    irccar = Irccar()
    irccar.read_text('2\n0.0 0.5\n1.0 0.5\n')
    assert(np.array_equal(irccar.points,[[0.0,0.5],[1.0,0.5]]))
    reread_irccar = Irccar()
    reread_irccar.read_text(irccar.write_text())
    assert(np.array_equal(reread_irccar.points,irccar.points))
#end def test_auxiliary_formatted_files


def test_extended_input_file_registry(tmp_path):
    from ..vasp_input import VaspInput,VRawFile

    raw_files = ('CHGCAR','DYNMATFULL','GAMMA','HESSEMAT','ML_AB','ML_FF',
                 'TAUCAR','WANPROJ')
    for filename in raw_files:
        (tmp_path / filename).write_text(filename+' synthetic data\n')
    #end for
    (tmp_path / 'KPOINTS_OPT').write_text(
        'mesh\n0\nGamma\n1 1 1\n'
        )

    vasp_input = VaspInput(tmp_path)
    for filename in raw_files:
        name = filename.lower()
        assert(isinstance(vasp_input[name],VRawFile))
        assert(vasp_input[name].text==filename+' synthetic data\n')
    #end for
    assert(vasp_input.kpoints_opt.mode=='auto')
    assert('EXHCAR' not in VaspInput.all_inputs)
#end def test_extended_input_file_registry


def test_potcar_construction_and_metadata(tmp_path):
    from ..vasp_input import Potcar

    empty = Potcar()
    assert(empty.write_text()=='')

    potcar_text = 'PAW_PBE H 01Jan2001\n1.0\nVRHFIN =H:\nEnd of Dataset'
    potcar_file = tmp_path / 'H.POTCAR'
    potcar_file.write_text(potcar_text)

    potcar = Potcar(potcar_file)
    assert(len(potcar.pseudos)==1)
    assert(potcar.write_text()==potcar_text)

    info = potcar.pot_info()
    assert(len(info)==1)
    assert(info[0].element=='H')
    assert(info[0].Zval==1)

    configured = Potcar(tmp_path,['H.POTCAR'])
    configured.load()
    assert(configured.pseudos[0]==potcar_text)
#end def test_potcar_construction_and_metadata


def test_input_filename_filtering(tmp_path):
    from ..vasp_input import VaspInput

    (tmp_path / 'INCAR').write_text('ENCUT = 100\n')
    (tmp_path / 'job_INCAR.in').write_text('ENCUT = 200\n')

    vasp_input = VaspInput(tmp_path,prefix='job_',postfix='.in')
    assert(vasp_input.incar.encut==200.0)
#end def test_input_filename_filtering



def test_empty_init():
    from ..vasp_input import VaspInput
    vi = VaspInput()
#end test_empty_init



def test_read():
    from ..vasp_input import VaspInput

    test_files_dir = TEST_FILES['diamond_INCAR'].parent
    vi = VaspInput(test_files_dir, prefix='diamond_')

    del vi.potcar.filepath

    check_vs_serial_reference(vi,'read')
#end test_read



def test_write(tmp_path):
    from ..vasp_input import VaspInput

    test_files_dir = TEST_FILES['diamond_INCAR'].parent

    vi_read = VaspInput(test_files_dir, prefix='diamond_')

    vi_read.write(tmp_path, prefix='write_diamond_')

    vi_write = VaspInput(tmp_path, prefix='write_diamond_')

    del vi_write.potcar.filepath

    check_vs_serial_reference(vi_write,'write')
#end test_write


@isolate_nexus_core
def test_generate(tmp_path):
    import numpy as np
    from ..nexus_base import nexus_noncore
    from ..physical_system import generate_physical_system
    from ..vasp_input import generate_vasp_input,VaspInput

    pseudo_dir = tmp_path / 'pseudopotentials'
    pseudo_dir.mkdir()

    nexus_core.local_directory  = str(tmp_path)
    nexus_core.remote_directory = str(tmp_path)
    nexus_core.file_locations = nexus_core.file_locations + [str(tmp_path)]
    nexus_noncore.pseudo_dir = pseudo_dir

    (pseudo_dir / 'C.POTCAR').write_text(c_potcar_text)
    PseudoSet.pseudo_files = {
        'C.POTCAR':str((pseudo_dir/'C.POTCAR').resolve())
        }


    dia16 = generate_physical_system(
        structure = TEST_FILES['d16bulk.POSCAR'],
        C         = 4
        )

    vi = generate_vasp_input(
        system       = dia16,
        pseudos      = ['C.POTCAR'],
        input_type   = 'generic',
        istart       = 0,
        icharg       = 2,
        encut        = 450,
        nsw          = 5,
        ibrion       = 2,
        isif         = 2,
        kcenter      = 'monkhorst',
        kgrid        = (2,2,2),
        kshift       = (0,0,0),
        )

    assert(isinstance(vi,VaspInput))

    del vi.potcar.filepath

    check_vs_serial_reference(vi,'generate')

    vi_minimal = generate_vasp_input(system=dia16)
    assert(vi_minimal.incar.nelect==64)

    vi_title = generate_vasp_input(title='diamond')
    assert(vi_title.incar.system=='diamond')

    kbasis = np.array([
        [0.25,0.00,0.00],
        [0.00,0.25,0.00],
        [0.00,0.00,0.25],
        ])
    kshift = (0.5,0.0,0.5)
    vi_basis = generate_vasp_input(
        kbasis=kbasis,
        kcoord='reciprocal',
        kshift=kshift,
        )
    assert(vi_basis.kpoints.mode=='basis')
    assert(vi_basis.kpoints.coord=='reciprocal')
    assert(np.array_equal(vi_basis.kpoints.kbasis,kbasis))
    assert(vi_basis.kpoints.kshift==kshift)
#end def test_generate


def test_return_system_sampling_and_electronic_state():
    import numpy as np
    from ..vasp_input import Incar,Kpoints,Poscar,Potcar,VaspInput

    poscar = Poscar()
    poscar.scale = 1.0
    poscar.axes = 2*np.eye(3)
    poscar.elem = None
    poscar.elem_count = np.array([1])
    poscar.coord = 'direct'
    poscar.pos = np.array([[0.0,0.0,0.0]])

    potcar = Potcar()
    potcar.read_text(
        'PAW_PBE H 01Jan2001\n'
        '1.0\n'
        'VRHFIN =H:\n'
        'End of Dataset\n'
        )
    kpoints = Kpoints()
    kpoints.mode = 'explicit'
    kpoints.coord = 'reciprocal'
    kpoints.kpoints = np.array([[0.25,0.0,0.0]])
    kpoints.kweights = np.array([2.0])

    vasp_input = VaspInput()
    vasp_input.poscar = poscar
    vasp_input.potcar = potcar
    vasp_input.kpoints = kpoints
    vasp_input.incar = Incar()
    vasp_input.incar.assign(nelect=0.5,nupdown=1)
    system = vasp_input.return_system()

    assert(np.array_equal(system.structure.elem,['H']))
    assert(np.allclose(system.structure.kpoints,[[np.pi/4,0,0]]))
    assert(np.array_equal(system.structure.kweights,[2.0]))
    assert(system.net_charge==0.5)
    assert(system.net_spin==1)
    assert(system.valency.H==1)

    del vasp_input.kpoints
    vasp_input.incar = Incar()
    vasp_input.incar.assign(kspacing=2.0,kgamma=False)
    system = vasp_input.return_system()
    assert(len(system.structure.kpoints)==8)

    vasp_input.kpoints = Kpoints()
    vasp_input.kpoints.read_text('automatic\n0\nAuto\n4.0\n')
    system = vasp_input.return_system()
    assert(len(system.structure.kpoints)==8)
#end def test_return_system_sampling_and_electronic_state


def test_generate_poscar_reorders_atom_data_and_accepts_system():
    import numpy as np
    from ..physical_system import PhysicalSystem
    from ..structure import Structure
    from ..vasp_input import generate_poscar,Poscar

    structure = Structure(
        axes   = np.eye(3),
        elem   = ['C','H','C'],
        pos    = [[0,0,0],[0.5,0,0],[0.25,0,0]],
        vel    = [[1,0,0],[2,0,0],[3,0,0]],
        frozen = [[True,False,False],
                  [False,True,False],
                  [False,False,True]],
        units  = 'A',
        )
    system = PhysicalSystem(structure=structure)

    poscar = generate_poscar(system,coord='direct')

    assert(poscar.elem==['C','H'])
    assert(poscar.elem_count==[2,1])
    assert(np.array_equal(poscar.pos[:,0],[0.0,0.25,0.5]))
    assert(poscar.vel_coord=='cartesian')
    assert(np.array_equal(poscar.vel[:,0],[1,3,2]))
    assert(np.array_equal(
        poscar.dynamic,
        [[False,True,True],
         [True,True,False],
         [True,False,True]],
        ))
    assert(np.array_equal(structure.elem,['C','H','C']))

    reread = Poscar()
    reread.read_text(poscar.write_text())
    assert(reread.vel_coord=='cartesian')
    assert(np.array_equal(reread.vel[:,0],[1,3,2]))
#end def test_generate_poscar_reorders_atom_data_and_accepts_system


def test_poscar_inline_comments_and_named_vector_block():
    import numpy as np
    from ..vasp_input import Poscar

    text = (
        'synthetic triclinic structure\n'
        '1.2345 ! arbitrary scale\n'
        '2.1 0.2 0.3 ! first lattice vector\n'
        '0.4 3.2 0.5 # second lattice vector\n'
        '0.6 0.7 4.3 ! third lattice vector\n'
        'Li Ne ! synthetic species pair\n'
        '2 1 # deliberately unequal counts\n'
        'Direct ! coordinate mode\n'
        '0.11 0.22 0.33 ! Li 1\n'
        '0.44 0.55 0.66 # Li 2\n'
        '0.77 0.88 0.99 ! Ne\n'
        'initial forces ! named vector block\n'
        '1.1 -2.2 3.3 ! vector 1\n'
        '-4.4 5.5 -6.6 ! vector 2\n'
        '7.7 -8.8 9.9 ! vector 3\n'
        )
    poscar = Poscar()
    poscar.read_text(text)

    assert(poscar.scale==1.2345)
    assert(np.allclose(
        poscar.axes,
        [[2.1,0.2,0.3],[0.4,3.2,0.5],[0.6,0.7,4.3]],
        ))
    assert(poscar.coord=='direct')
    assert(poscar.vel is None)
    assert(poscar.vector_header=='initial forces')
    assert(np.allclose(
        poscar.vectors,
        [[1.1,-2.2,3.3],[-4.4,5.5,-6.6],[7.7,-8.8,9.9]],
        ))

    reread = Poscar()
    reread.read_text(poscar.write_text())
    assert(reread.vector_header=='initial forces')
    assert(np.array_equal(reread.vectors,poscar.vectors))
#end def test_poscar_inline_comments_and_named_vector_block


def test_vasp_input_run_type():
    from ..vasp_input import Incar,VaspInput

    vasp_input = VaspInput()
    assert(vasp_input.run_type()=='unknown')

    cases = (
        (dict(),                         'static'),
        (dict(nsw=10),                   'static'),
        (dict(ibrion=0,nsw=10),          'md'),
        (dict(ibrion=2,nsw=10),          'relax'),
        (dict(ibrion=6,nsw=1),           'phonon'),
        (dict(ibrion=2,nsw=10,images=3), 'neb'),
    )
    for inputs,expected in cases:
        vasp_input.incar = Incar()
        vasp_input.incar.assign(**inputs)
        assert(vasp_input.run_type()==expected
               )
    #end for
#end def test_vasp_input_run_type


def test_neb_input_roundtrip(tmp_path):
    import numpy as np
    from ..vasp_input import Incar,NebPoscars,Poscar,VaspInput

    vasp_input = VaspInput()
    vasp_input.incar = Incar()
    vasp_input.incar.images = 1
    vasp_input.poscar = NebPoscars()
    for n,x in enumerate((0.0,0.5,1.0)):
        poscar = Poscar()
        poscar.scale = 1.0
        poscar.axes = np.eye(3)
        poscar.elem = np.array(['H'])
        poscar.elem_count = np.array([1])
        poscar.coord = 'direct'
        poscar.pos = np.array([[x,0.0,0.0]])
        vasp_input.poscar[n] = poscar
    #end for

    vasp_input.write(tmp_path)
    reread = VaspInput(tmp_path)

    assert(isinstance(reread.poscar,NebPoscars))
    assert(sorted(reread.poscar.keys())==[0,1,2])
    assert(np.allclose(reread.poscar[1].pos,[[0.5,0.0,0.0]]))
#end def test_neb_input_roundtrip


def test_vasp_input_validation():
    import numpy as np
    from ..vasp_input import Incar,Kpoints,Poscar,Potcar,VaspInput

    vasp_input = VaspInput()
    vasp_input.incar = Incar()
    vasp_input.poscar = Poscar()
    vasp_input.poscar.scale = 1.0
    vasp_input.poscar.axes = np.eye(3)
    vasp_input.poscar.elem = np.array(['H'])
    vasp_input.poscar.elem_count = np.array([1])
    vasp_input.poscar.coord = 'direct'
    vasp_input.poscar.pos = np.array([[0.0,0.0,0.0]])
    vasp_input.potcar = Potcar()
    vasp_input.potcar.read_text(
        'PAW_PBE H 01Jan2001\n'
        '1.0\n'
        'VRHFIN =H:\n'
        'End of Dataset\n'
        )

    assert(vasp_input.validate(exit=False))

    vasp_input.kpoints = Kpoints()
    vasp_input.kpoints.mode = 'explicit'
    vasp_input.kpoints.coord = 'reciprocal'
    vasp_input.kpoints.kpoints = np.zeros((2,3))
    vasp_input.kpoints.kweights = np.ones(1)
    assert(not vasp_input.validate(exit=False))
    with pytest.raises(NexusError,match='KPOINTS weights'):
        vasp_input.validate()
    #end with
#end def test_vasp_input_validation
