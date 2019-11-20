
import testing
from testing import value_eq,object_eq,text_eq


def test_import():
    from gamess_analyzer import GamessAnalyzer
#end def test_import


def test_empty_init():
    from gamess_analyzer import GamessAnalyzer

    GamessAnalyzer()
#end def test_empty_init



def test_analyze():
    import os
    from numpy import array,ndarray
    from generic import obj
    from gamess_analyzer import GamessAnalyzer

    tpath = testing.setup_unit_test_output_directory(
        test      = 'gamess_analyzer',
        subtest   = 'test_analyze',
        file_sets = ['gms.inp','gms.out'],
        )

    outpath = os.path.join(tpath,'gms.inp')

    # no analysis
    ga = GamessAnalyzer(outpath)

    del ga.info.input
    del ga.info.path

    ga_ref = obj(
        info = obj(
            exit            = False,
            initialized     = True,
            prefix          = 'gms',
            files = obj(
                aabb41          = 'gms.F41',
                aoints          = 'gms.F08',
                bbaa42          = 'gms.F42',
                bbbb43          = 'gms.F43',
                casints         = 'gms.F13',
                ccdiis          = 'gms.F71',
                ccints          = 'gms.F72',
                ccquads         = 'gms.F78',
                ccrest          = 'gms.F70',
                cct1amp         = 'gms.F73',
                cct2amp         = 'gms.F74',
                cct3amp         = 'gms.F75',
                ccve            = 'gms.F77',
                ccvm            = 'gms.F76',
                ciints          = 'gms.F14',
                civectr         = 'gms.F12',
                csfsave         = 'gms.F17',
                dafl30          = 'gms.F30',
                dasort          = 'gms.F20',
                dcphfh2         = 'gms.F67',
                dftgrid         = 'gms.F22',
                dftints         = 'gms.F21',
                dictnry         = 'gms.F10',
                drtfile         = 'gms.F11',
                efmof           = 'gms.F103',
                efmoi           = 'gms.F102',
                efpind          = 'gms.F25',
                eomdg12         = 'gms.F89',
                eomhc1          = 'gms.F83',
                eomhc2          = 'gms.F84',
                eomhhhh         = 'gms.F85',
                eomhl1          = 'gms.F98',
                eomhl2          = 'gms.F99',
                eomlvec         = 'gms.F97',
                eompppp         = 'gms.F86',
                eomramp         = 'gms.F87',
                eomrtmp         = 'gms.F88',
                eomstar         = 'gms.F80',
                eomvec1         = 'gms.F81',
                eomvec2         = 'gms.F82',
                eomvl1          = 'gms.F95',
                eomvl2          = 'gms.F96',
                fockder         = 'gms.F18',
                hessian         = 'gms.F38',
                input           = 'gms.inp',
                jkfile          = 'gms.F23',
                mcqd50          = 'gms.F50',
                mcqd51          = 'gms.F51',
                mcqd52          = 'gms.F52',
                mcqd53          = 'gms.F53',
                mcqd54          = 'gms.F54',
                mcqd55          = 'gms.F55',
                mcqd56          = 'gms.F56',
                mcqd57          = 'gms.F57',
                mcqd58          = 'gms.F58',
                mcqd59          = 'gms.F59',
                mcqd60          = 'gms.F60',
                mcqd61          = 'gms.F61',
                mcqd62          = 'gms.F62',
                mcqd63          = 'gms.F63',
                mcqd64          = 'gms.F64',
                mltpl           = 'gms.F28',
                mltplt          = 'gms.F29',
                mmciitr         = 'gms.F94',
                mmcivc1         = 'gms.F93',
                mmcivec         = 'gms.F92',
                mmhpp           = 'gms.F91',
                mmpp            = 'gms.F90',
                moints          = 'gms.F09',
                nmrint1         = 'gms.F61',
                ordint          = 'gms.F24',
                output          = 'gms.out',
                pcmdata         = 'gms.F26',
                pcmints         = 'gms.F27',
                punch           = 'gms.F07',
                quadsvo         = 'gms.F79',
                remd            = 'gms.F44',
                restart         = 'gms.F35',
                soccdat         = 'gms.F40',
                traject         = 'gms.F04',
                work15          = 'gms.F15',
                work16          = 'gms.F16',
                work19          = 'gms.F19',
                ),
            )
        )

    assert(object_eq(ga.to_obj(),ga_ref))


    # full analysis
    ga = GamessAnalyzer(outpath,analyze=True)

    del ga.info.input
    del ga.info.path
    del ga.info.files

    ga_ref = obj(
        ao_populations = obj(
            lowdin          = array(
                              [8.3070e-02, 2.5080e-02, 5.2000e-04, 2.3000e-04, 7.3660e-02, 0.0000e+00,
                               0.0000e+00, 1.3479e-01, 0.0000e+00, 0.0000e+00, 6.2610e-02, 0.0000e+00,
                               0.0000e+00, 2.1910e-02, 3.4700e-02, 3.4700e-02, 1.2332e-01, 0.0000e+00,
                               0.0000e+00, 0.0000e+00, 2.3050e-02, 2.3050e-02, 1.2369e-01, 0.0000e+00,
                               0.0000e+00, 0.0000e+00, 0.0000e+00, 0.0000e+00, 9.8720e-02, 0.0000e+00,
                               1.6890e-02, 0.0000e+00, 1.6890e-02, 0.0000e+00, 0.0000e+00, 0.0000e+00,
                               3.0345e-01, 1.7130e-02, 7.4000e-04, 7.0500e-03, 2.4783e-01, 0.0000e+00,
                               0.0000e+00, 5.0000e-05, 0.0000e+00, 0.0000e+00, 1.1000e-03, 0.0000e+00,
                               0.0000e+00, 4.7200e-03, 6.1080e-02, 6.1080e-02, 6.1380e-02, 0.0000e+00,
                               0.0000e+00, 0.0000e+00, 1.0230e-01, 1.0230e-01, 1.3176e-01, 0.0000e+00,
                               0.0000e+00, 0.0000e+00, 0.0000e+00, 0.0000e+00, 7.6000e-04, 0.0000e+00,
                               2.1000e-04, 0.0000e+00, 2.1000e-04, 0.0000e+00, 0.0000e+00, 0.0000e+00],
                               dtype=float),
            lowdin_shell    = array(
                              [8.3070e-02, 2.5080e-02, 5.2000e-04, 2.3000e-04, 7.3660e-02, 3.0345e-01,
                               1.7130e-02, 7.4000e-04, 7.0500e-03, 2.4783e-01, 1.3479e-01, 6.2610e-02,
                               2.1910e-02, 5.0000e-05, 1.1000e-03, 4.7200e-03, 1.9272e-01, 1.6979e-01,
                               1.8354e-01, 3.3636e-01, 1.3250e-01, 1.1800e-03],
                               dtype=float),
            mulliken        = array(
                              [ 7.4225e-01,  1.1141e-01, -1.0500e-03, -4.9000e-04, -4.5679e-01,  0.0000e+00,
                                0.0000e+00,  4.6670e-01,  0.0000e+00,  0.0000e+00, -2.0338e-01,  0.0000e+00,
                                0.0000e+00,  2.5550e-02,  2.3860e-02,  0.0000e+00,  0.0000e+00,  0.0000e+00,
                                0.0000e+00,  0.0000e+00,  1.4390e-02,  0.0000e+00,  0.0000e+00,  0.0000e+00,
                                0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00,  6.2800e-03,
                                0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00,
                                1.5046e+00,  2.9270e-02,  1.8600e-03,  2.0000e-05, -2.6496e-01,  0.0000e+00,
                                0.0000e+00,  3.0000e-05,  0.0000e+00,  0.0000e+00,  8.3000e-04,  0.0000e+00,
                                0.0000e+00, -4.2000e-04, -0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00,
                                0.0000e+00,  0.0000e+00,  6.0000e-05,  0.0000e+00,  0.0000e+00,  0.0000e+00,
                                0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00, -0.0000e+00,
                                0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00,  0.0000e+00],
                                dtype=float),
            mulliken_shell  = array(
                              [ 7.4225e-01,  1.1141e-01, -1.0500e-03, -4.9000e-04, -4.5679e-01,  1.5046e+00,
                                2.9270e-02,  1.8600e-03,  2.0000e-05, -2.6496e-01,  4.6670e-01, -2.0338e-01,
                                2.5550e-02,  3.0000e-05,  8.3000e-04, -4.2000e-04,  2.3860e-02,  1.4390e-02,
                                0.0000e+00,  6.0000e-05,  6.2800e-03,  0.0000e+00],
                                dtype=float),
            basis = obj(
                angular         = array(
                                   ['S','S','S','S','S','X','Y','Z','X','Y','Z','X','Y','Z','XX','YY','ZZ',
                                   'XY','XZ','YZ','XX','YY','ZZ','XY','XZ','YZ','XXX','YYY','ZZZ','XXY',
                                   'XXZ','YYX','YYZ','ZZX','ZZY','XYZ','S','S','S','S','S','X','Y','Z','X',
                                   'Y','Z','X','Y','Z','XX','YY','ZZ','XY','XZ','YZ','XX','YY','ZZ','XY',
                                   'XZ','YZ','XXX','YYY','ZZZ','XXY','XXZ','YYX','YYZ','ZZX','ZZY','XYZ']),
                d               = array([14,15,16,17,18,19,20,21,22,23,24,25,50,51,52,53,54,55,56,57,58,59,60,61],dtype=int),
                element         = array(
                                  ['Li','Li','Li','Li','Li','Li','Li','Li','Li','Li','Li','Li','Li','Li',
                                   'Li','Li','Li','Li','Li','Li','Li','Li','Li','Li','Li','Li','Li','Li',
                                   'Li','Li','Li','Li','Li','Li','Li','Li','H','H','H','H','H','H','H','H',
                                   'H','H','H','H','H','H','H','H','H','H','H','H','H','H','H','H','H','H',
                                   'H','H','H','H','H','H','H','H','H','H']),
                f               = array([26,27,28,29,30,31,32,33,34,35,62,63,64,65,66,67,68,69,70,71],dtype=int),
                g               = array([],dtype=int),
                p               = array([ 5,6,7,8,9,10,11,12,13,41,42,43,44,45,46,47,48,49],dtype=int),
                s               = array([ 0,1,2,3,4,36,37,38,39,40],dtype=int),
                spec_index      = array([1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,2,
                                   2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2],dtype=int),
                ),
            lowdin_angular = obj(
                d               = 0.88241,
                f               = 0.13368,
                g               = 0.0,
                p               = 0.22518,
                s               = 0.75876,
                ),
            mulliken_angular = obj(
                d               = 0.03831,
                f               = 0.00628,
                g               = 0.0,
                p               = 0.28931,
                s               = 1.666119999999999,
                ),
            shell = obj(
                d               = array([16,17,18,19],dtype=int),
                f               = array([20,21],dtype=int),
                g               = array([],dtype=int),
                p               = array([10,11,12,13,14,15],dtype=int),
                s               = array([0,1,2,3,4,5,6,7,8,9],dtype=int),
                ),
            ),
        counts = obj(
            cao             = 72,
            mos             = 62,
            shells          = 22,
            ),
        energy = obj(
            electron_electron_potential = 0.4814819907,
            nuclear_repulsion           = 0.3317933721,
            nucleus_electron_potential  = -2.199450924,
            nucleus_nucleus_potential   = 0.3317933721,
            one_electron                = -1.5640017754,
            total                       = -0.7507264125,
            total_kinetic               = 0.6354491487,
            total_potential             = -1.3861755612,
            two_electron                = 0.4814819907,
            ),
        info = obj(
            exit            = False,
            initialized     = True,
            prefix          = 'gms',
            ),
        )

    assert(object_eq(ga.to_obj(),ga_ref))

#end def test_analyze
