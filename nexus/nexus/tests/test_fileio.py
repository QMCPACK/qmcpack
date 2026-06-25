import pytest
from . import NexusTestOrder
pytestmark = pytest.mark.order(NexusTestOrder.FILEIO)

from ..generic import generic_settings
generic_settings.raise_error = True

from . import TEST_DIR
from ..testing import value_eq, object_eq


TEST_FILES = {
    "scf.in":               TEST_DIR / "test_fileio_files/scf.in",
    "VO2_R_48_dens.CHGCAR": TEST_DIR / "test_fileio_files/VO2_R_48_dens.CHGCAR",
    "VO2_R_48_dens.xsf":    TEST_DIR / "test_fileio_files/VO2_R_48_dens.xsf",
    "VO2_R_48.POSCAR":      TEST_DIR / "test_fileio_files/VO2_R_48.POSCAR",
    "VO2_R_48.xsf":         TEST_DIR / "test_fileio_files/VO2_R_48.xsf",
    }

for file in TEST_FILES.values():
    assert(file.exists()), f"Test file not found! {file}"


def test_textfile():
    from ..fileio import TextFile

    # test empty initialization
    empty = TextFile()

    # test read
    f = TextFile(TEST_FILES["scf.in"])

    assert(len(f.read())==1225)
    assert(len(f.lines())==55)
    assert(len(f.tokens())==125)

    assert(f.readline()=='&CONTROL\n')
    assert(f.readline()=="   calculation     = 'scf'\n")

    assert(f.readline('celldm')=='celldm(1)       = 1.0\n')

    assert(f.readtokens('ecutwfc')==['ecutwfc', '=', '272'])

    val = f.readtokensf('CELL_PARAMETERS',float,float,float)
    ref = [28.34589199, 0.0, 0.0]
    assert(value_eq(val,ref))

    f.close()
#end def test_textfile


def test_xsffile(tmp_path):
    import numpy as np
    from ..fileio import XsfFile

    # test empty initialization
    empty = XsfFile()
    assert(not empty.is_valid())

    # populate reference object
    ref = XsfFile()
    ref.set(
        elem = np.array([8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,
                         8,8,8,8,8,8,8,8,8,8,8,8,23,23,23,23,23,
                         23,23,23,23,23,23,23,23,23,23,23],dtype=int),
        filetype    = 'xsf',
        periodicity = 'crystal',
        pos         = np.array([
                [ 1.36683546,  1.36683546,  0.        ],
                [ 3.18776454,  3.18776454,  0.        ],
                [ 3.64413546,  0.91046454,  1.4264    ],
                [ 5.46506454, -0.91046454,  1.4264    ],
                [ 1.36683546,  1.36683546,  2.8528    ],
                [ 3.18776454,  3.18776454,  2.8528    ],
                [ 3.64413546,  0.91046454,  4.2792    ],
                [ 5.46506454, -0.91046454,  4.2792    ],
                [ 5.92143546,  1.36683546,  0.        ],
                [ 3.18776454, -1.36683546,  0.        ],
                [ 3.64413546, -3.64413546,  1.4264    ],
                [ 0.91046454, -0.91046454,  1.4264    ],
                [ 5.92143546,  1.36683546,  2.8528    ],
                [ 3.18776454, -1.36683546,  2.8528    ],
                [ 3.64413546, -3.64413546,  4.2792    ],
                [ 0.91046454, -0.91046454,  4.2792    ],
                [ 1.36683546,  1.36683546,  5.7056    ],
                [ 3.18776454,  3.18776454,  5.7056    ],
                [ 3.64413546,  0.91046454,  7.132     ],
                [ 5.46506454, -0.91046454,  7.132     ],
                [ 1.36683546,  1.36683546,  8.5584    ],
                [ 3.18776454,  3.18776454,  8.5584    ],
                [ 3.64413546,  0.91046454,  9.9848    ],
                [ 5.46506454, -0.91046454,  9.9848    ],
                [ 5.92143546,  1.36683546,  5.7056    ],
                [ 3.18776454, -1.36683546,  5.7056    ],
                [ 3.64413546, -3.64413546,  7.132     ],
                [ 0.91046454, -0.91046454,  7.132     ],
                [ 5.92143546,  1.36683546,  8.5584    ],
                [ 3.18776454, -1.36683546,  8.5584    ],
                [ 3.64413546, -3.64413546,  9.9848    ],
                [ 0.91046454, -0.91046454,  9.9848    ],
                [ 0.        ,  0.        ,  0.        ],
                [ 2.2773    ,  2.2773    ,  1.4264    ],
                [ 0.        ,  0.        ,  2.8528    ],
                [ 2.2773    ,  2.2773    ,  4.2792    ],
                [ 4.5546    ,  0.        ,  0.        ],
                [ 2.2773    , -2.2773    ,  1.4264    ],
                [ 4.5546    ,  0.        ,  2.8528    ],
                [ 2.2773    , -2.2773    ,  4.2792    ],
                [ 0.        ,  0.        ,  5.7056    ],
                [ 2.2773    ,  2.2773    ,  7.132     ],
                [ 0.        ,  0.        ,  8.5584    ],
                [ 2.2773    ,  2.2773    ,  9.9848    ],
                [ 4.5546    ,  0.        ,  5.7056    ],
                [ 2.2773    , -2.2773    ,  7.132     ],
                [ 4.5546    ,  0.        ,  8.5584    ],
                [ 2.2773    , -2.2773    ,  9.9848    ]
                ],dtype=float),
        primvec = np.array([
                [  4.5546,  -4.5546,   0.    ],
                [  4.5546,   4.5546,   0.    ],
                [  0.    ,   0.    ,  11.4112]
                ],dtype=float),
        )
    assert(ref.is_valid())

    # test read
    f = XsfFile(TEST_FILES['VO2_R_48.xsf'])
    assert(f.is_valid())
    assert(object_eq(f,ref))

    # test write
    outfile = tmp_path / "test.xsf"
    f.write(outfile)
    f2 = XsfFile(outfile)
    assert(f2.is_valid())
    assert(object_eq(f2,ref))
#end def test_xsffile


def test_xsffile_density(tmp_path):
    import numpy as np
    from ..fileio import XsfFile

    ref = XsfFile(TEST_FILES['VO2_R_48.xsf'])

    grid = 3,5,7
    dens = 0.01*np.arange(np.prod(grid),dtype=float)
    dens.shape=grid

    ref.add_density(ref.primvec,dens,add_ghost=True)
    assert(ref.is_valid())

    f = XsfFile(TEST_FILES['VO2_R_48_dens.xsf'])
    assert(f.is_valid())
    assert(object_eq(f,ref))

    d = f.get_density().values
    assert(isinstance(d,np.ndarray))
    assert(d.shape==(4,6,8))

    outfile = tmp_path / "test_density.xsf"
    f.write(outfile)
    f2 = XsfFile(outfile)
    assert(f2.is_valid())
    assert(object_eq(f2,ref))
#end def test_xsffile_density



def test_poscar_file(tmp_path):
    import numpy as np
    from ..fileio import PoscarFile

    # test empty initialization
    empty = PoscarFile()
    assert(not empty.is_valid())

    # populate reference object
    ref = PoscarFile()
    ref.set(
        axes        = np.array([
                [  4.5546,  -4.5546,   0.    ],
                [  4.5546,   4.5546,   0.    ],
                [  0.    ,   0.    ,  11.4112]
                ],dtype=float),
        coord       = 'cartesian',
        description = 'System cell and coordinates',
        dynamic     = None,
        elem        = np.array(['O', 'V'],dtype=str),
        elem_count  = np.array([32, 16],dtype=int),
        pos         = np.array([
                [ 1.36683546,  1.36683546,  0.        ],
                [ 3.18776454,  3.18776454,  0.        ],
                [ 3.64413546,  0.91046454,  1.4264    ],
                [ 5.46506454, -0.91046454,  1.4264    ],
                [ 1.36683546,  1.36683546,  2.8528    ],
                [ 3.18776454,  3.18776454,  2.8528    ],
                [ 3.64413546,  0.91046454,  4.2792    ],
                [ 5.46506454, -0.91046454,  4.2792    ],
                [ 5.92143546,  1.36683546,  0.        ],
                [ 3.18776454, -1.36683546,  0.        ],
                [ 3.64413546, -3.64413546,  1.4264    ],
                [ 0.91046454, -0.91046454,  1.4264    ],
                [ 5.92143546,  1.36683546,  2.8528    ],
                [ 3.18776454, -1.36683546,  2.8528    ],
                [ 3.64413546, -3.64413546,  4.2792    ],
                [ 0.91046454, -0.91046454,  4.2792    ],
                [ 1.36683546,  1.36683546,  5.7056    ],
                [ 3.18776454,  3.18776454,  5.7056    ],
                [ 3.64413546,  0.91046454,  7.132     ],
                [ 5.46506454, -0.91046454,  7.132     ],
                [ 1.36683546,  1.36683546,  8.5584    ],
                [ 3.18776454,  3.18776454,  8.5584    ],
                [ 3.64413546,  0.91046454,  9.9848    ],
                [ 5.46506454, -0.91046454,  9.9848    ],
                [ 5.92143546,  1.36683546,  5.7056    ],
                [ 3.18776454, -1.36683546,  5.7056    ],
                [ 3.64413546, -3.64413546,  7.132     ],
                [ 0.91046454, -0.91046454,  7.132     ],
                [ 5.92143546,  1.36683546,  8.5584    ],
                [ 3.18776454, -1.36683546,  8.5584    ],
                [ 3.64413546, -3.64413546,  9.9848    ],
                [ 0.91046454, -0.91046454,  9.9848    ],
                [ 0.        ,  0.        ,  0.        ],
                [ 2.2773    ,  2.2773    ,  1.4264    ],
                [ 0.        ,  0.        ,  2.8528    ],
                [ 2.2773    ,  2.2773    ,  4.2792    ],
                [ 4.5546    ,  0.        ,  0.        ],
                [ 2.2773    , -2.2773    ,  1.4264    ],
                [ 4.5546    ,  0.        ,  2.8528    ],
                [ 2.2773    , -2.2773    ,  4.2792    ],
                [ 0.        ,  0.        ,  5.7056    ],
                [ 2.2773    ,  2.2773    ,  7.132     ],
                [ 0.        ,  0.        ,  8.5584    ],
                [ 2.2773    ,  2.2773    ,  9.9848    ],
                [ 4.5546    ,  0.        ,  5.7056    ],
                [ 2.2773    , -2.2773    ,  7.132     ],
                [ 4.5546    ,  0.        ,  8.5584    ],
                [ 2.2773    , -2.2773    ,  9.9848    ]
                ],dtype=float),
        scale       = 1.0,
        vel         = None,
        vel_coord   = None,
        )

    # test read
    f = PoscarFile(TEST_FILES['VO2_R_48.POSCAR'])
    assert(f.is_valid())
    assert(object_eq(f,ref))

    # test write
    outfile = tmp_path / "test.POSCAR"
    f.write(outfile)
    f2 = PoscarFile(outfile)
    assert(f2.is_valid())
    assert(object_eq(f2,ref))

    # test incorporate xsf
    from ..fileio import XsfFile
    x = XsfFile(TEST_FILES['VO2_R_48.xsf'])
    f = PoscarFile()
    f.incorporate_xsf(x)
    assert(f.is_valid())
    assert(object_eq(f,ref))
#end def test_poscar_file



def test_chgcar_file(tmp_path):
    from ..fileio import XsfFile
    from ..fileio import PoscarFile
    from ..fileio import ChgcarFile

    empty = ChgcarFile()
    assert(not empty.is_valid())

    # get reference poscar and xsf files
    p = PoscarFile(TEST_FILES['VO2_R_48.POSCAR'])
    x = XsfFile(TEST_FILES['VO2_R_48_dens.xsf'])

    # create and test reference chgcar file
    ref = ChgcarFile()
    ref.incorporate_xsf(x)
    assert(ref.is_valid())
    assert(object_eq(ref.poscar,p))

    # test read
    f = ChgcarFile(TEST_FILES['VO2_R_48_dens.CHGCAR'])
    assert(f.is_valid())
    assert(object_eq(f,ref))

    # test write
    outfile = tmp_path / "test.CHGCAR"
    f.write(outfile)
    f2 = ChgcarFile(outfile)
    assert(f2.is_valid())
    assert(object_eq(f2,ref))
#end def test_chgcar_file
