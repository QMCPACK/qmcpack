
import testing
from testing import value_eq,object_eq


associated_files = dict()

def get_files():
    return testing.collect_unit_test_file_paths('vasp_analyzer',associated_files)
#end def get_files


relax_files = [
    'relax.CONTCAR',
    'relax.DOSCAR',
    'relax.EIGENVAL',
    'relax.err',
    'relax.IBZKPT',
    'relax.INCAR',
    'relax.KPOINTS',
    'relax.OSZICAR',
    'relax.out',
    'relax.OUTCAR',
    'relax.PCDAT',
    'relax.POSCAR',
    'relax.qsub.in',
    'relax.vasprun.xml',
    'relax.XDATCAR',
    ]


def test_files():
    filenames = relax_files
    files = get_files()
    assert(set(filenames)==set(files.keys()))
#end def test_files



def test_import():
    from vasp_analyzer import VaspAnalyzer,VXML,OutcarData
#end def test_import



def test_empty_init():
    from vasp_analyzer import VaspAnalyzer

    va = VaspAnalyzer(None)
#end def test_empty_init



def test_analyze():
    import os
    from numpy import array,ndarray
    from generic import obj
    from vasp_analyzer import VaspAnalyzer,VXML,VXMLcoll,OutcarData

    tpath = testing.setup_unit_test_output_directory(
        test      = 'vasp_analyzer',
        subtest   = 'test_analyze',
        file_sets = relax_files,
        )

    incar_path = os.path.join(tpath,'relax.INCAR')


    # empty init
    va = VaspAnalyzer(incar_path)

    va_ref = obj(
        info = obj(
            incar_file      = 'relax.INCAR',
            neb             = False,
            outcar_file     = 'relax.OUTCAR',
            prefix          = 'relax.',
            xml             = False,
            xml_file        = 'relax.vasprun.xml',
            incar = obj(
                ediff           = 1e-06,
                ediffg          = -0.01,
                elmin           = 5,
                ibrion          = 1,
                istart          = 0,
                nelect          = 130.0,
                nsw             = 20,
                prec            = 'normal',
                ),
            ),
        )

    del va.info.path

    assert(object_eq(va,va_ref))


    # full analysis
    va = VaspAnalyzer(incar_path,xml=True,analyze=True)

    del va.info.path

    info = obj(
        incar_file      = 'relax.INCAR',
        neb             = False,
        outcar_file     = 'relax.OUTCAR',
        prefix          = 'relax.',
        xml             = True,
        xml_file        = 'relax.vasprun.xml',
        incar = obj(
            ediff           = 1e-06,
            ediffg          = -0.01,
            elmin           = 5,
            ibrion          = 1,
            istart          = 0,
            nelect          = 130.0,
            nsw             = 20,
            prec            = 'normal',
            ),
        )
    assert(object_eq(va.info,info))

    types = dict(
        Efermi               = float  ,               
        core_potential_radii = ndarray,
        core_potentials      = ndarray,
        force                = ndarray,
        info                 = obj    ,
        ion_steps            = obj    ,
        lattice_vectors      = ndarray,
        memory               = obj    ,
        position             = ndarray,
        pressure             = float  ,
        stress               = obj    ,
        stress_kb            = obj    ,
        time                 = obj    ,
        total_drift          = ndarray,
        total_energy         = float  ,
        volume               = float  ,
        xmldata              = VXML   ,
        )
    for name,type in types.items():
        assert(name in va)
        assert(isinstance(va[name],type))
    #end for

    simple_data = obj(
        Efermi               = -0.5064,
        core_potential_radii = array([1.2059],dtype=float),
        core_potentials = array([-56.6601,-56.6601,-56.6601,-56.6601,-57.1969,
                                 -57.2036,-57.5372,-57.1969,-56.9801,-56.9801,
                                 -56.9801,-56.9801,-56.4963],dtype=float),
        force           = array(
                          [[ 0.132123, -0.132123, -0.009177],
                           [-0.132123, -0.132123, -0.009177],
                           [ 0.132123,  0.132123, -0.009177],
                           [-0.132123,  0.132123, -0.009177],
                           [ 0.      ,  0.      ,  0.002104],
                           [ 0.      ,  0.      , -0.000878],
                           [ 0.      ,  0.      ,  0.004689],
                           [ 0.      ,  0.      ,  0.002104],
                           [-0.0019  ,  0.0019  ,  0.005077],
                           [ 0.0019  ,  0.0019  ,  0.005077],
                           [-0.0019  , -0.0019  ,  0.005077],
                           [ 0.0019  , -0.0019  ,  0.005077],
                           [ 0.      ,  0.      ,  0.008385]],dtype=float),
        lattice_vectors = array(
                          [[ 5.62024,  0.     ,  0.     ],
                           [ 0.     ,  5.62024,  0.     ],
                           [ 0.     ,  0.     , 16.86072]],dtype=float),
        position        = array(
                          [[1.40506, 1.40506, 1.98704],
                           [4.21518, 1.40506, 1.98704],
                           [1.40506, 4.21518, 1.98704],
                           [4.21518, 4.21518, 1.98704],
                           [0.     , 0.     , 3.88944],
                           [0.     , 2.81012, 3.89432],
                           [2.81012, 0.     , 3.95815],
                           [2.81012, 2.81012, 3.88944],
                           [1.44092, 1.3692 , 5.85353],
                           [4.17932, 1.3692 , 5.85353],
                           [1.44092, 4.25104, 5.85353],
                           [4.17932, 4.25104, 5.85353],
                           [0.     , 2.81012, 7.49145]],dtype=float),
        pressure        = -29.47,
        total_drift     = array([0.,       0.,       0.042555],dtype=float),
        total_energy    = -71.46392704,
        volume          = 532.58,
        )
    for name,data in simple_data.items():
        assert(value_eq(va[name],data))
    #end for

    memory = obj(
        average = 0.0,
        maximum = 58808.0,
        )
    assert(object_eq(va.memory,memory))

    stress = obj(
        xx = -12.73119,
        xy = 0.0,
        yy = -12.73119,
        yz = 0.0,
        zx = 0.0,
        zz = -3.92272,
        )
    assert(object_eq(va.stress,stress))

    stress_kb = obj(
        xx = -38.29957,
        xy = 0.0,
        yy = -38.29957,
        yz = 0.0,
        zx = 0.0,
        zz = -11.80082,
        )
    assert(object_eq(va.stress_kb,stress_kb))

    time = obj(
        cpu     = 229.665,
        elapsed = 245.015,
        system  = 5.815,
        user    = 223.85,
        )
    assert(object_eq(va.time,time))

    assert(len(va.ion_steps)==10)
    for n in range(1,11):
        assert(n in va.ion_steps)
    #end for

    last_step = va.ion_steps[10]
    assert(len(last_step)==4)
    for n in range(1,4):
        assert(n in last_step)
        assert(isinstance(last_step[n],OutcarData))
    #end for

    last_data = obj(
        Efermi          = -0.5064,
        core_potential_radii = array([1.2059],dtype=float),
        core_potentials = array([-56.6601,-56.6601,-56.6601,-56.6601,-57.1969,
                                 -57.2036,-57.5372,-57.1969,-56.9801,-56.9801,
                                 -56.9801,-56.9801,-56.4963],dtype=float),
        force           = array(
                      [[ 0.132123, -0.132123, -0.009177],
                       [-0.132123, -0.132123, -0.009177],
                       [ 0.132123,  0.132123, -0.009177],
                       [-0.132123,  0.132123, -0.009177],
                       [ 0.      ,  0.      ,  0.002104],
                       [ 0.      ,  0.      , -0.000878],
                       [ 0.      ,  0.      ,  0.004689],
                       [ 0.      ,  0.      ,  0.002104],
                       [-0.0019  ,  0.0019  ,  0.005077],
                       [ 0.0019  ,  0.0019  ,  0.005077],
                       [-0.0019  , -0.0019  ,  0.005077],
                       [ 0.0019  , -0.0019  ,  0.005077],
                       [ 0.      ,  0.      ,  0.008385]],dtype=float),
        lattice_vectors = array(
                          [[ 5.62024,  0.     ,  0.     ],
                           [ 0.     ,  5.62024,  0.     ],
                           [ 0.     ,  0.     , 16.86072]],dtype=float),
        position        = array(
                          [[1.40506, 1.40506, 1.98704],
                           [4.21518, 1.40506, 1.98704],
                           [1.40506, 4.21518, 1.98704],
                           [4.21518, 4.21518, 1.98704],
                           [0.     , 0.     , 3.88944],
                           [0.     , 2.81012, 3.89432],
                           [2.81012, 0.     , 3.95815],
                           [2.81012, 2.81012, 3.88944],
                           [1.44092, 1.3692 , 5.85353],
                           [4.17932, 1.3692 , 5.85353],
                           [1.44092, 4.25104, 5.85353],
                           [4.17932, 4.25104, 5.85353],
                           [0.     , 2.81012, 7.49145]],dtype=float),
        pressure        = -29.47,
        total_drift     = array([0.,       0.,       0.042555],dtype=float),
        total_energy    = -71.46392704,
        volume          = 532.58,
        memory = obj(
            average         = 0.0,
            maximum         = 58808.0,
            ),
        stress = obj(
            xx              = -12.73119,
            xy              = 0.0,
            yy              = -12.73119,
            yz              = 0.0,
            zx              = 0.0,
            zz              = -3.92272,
            ),
        stress_kb = obj(
            xx              = -38.29957,
            xy              = 0.0,
            yy              = -38.29957,
            yz              = 0.0,
            zx              = 0.0,
            zz              = -11.80082,
            ),
        time = obj(
            cpu             = 229.665,
            elapsed         = 245.015,
            system          = 5.815,
            user            = 223.85,
            ),
        )

    assert(object_eq(va.ion_steps[10][4],last_data))
    

    vxml = va.xmldata

    keys = '''
        atominfo
        calculation
        finalpos
        generator
        incar
        initialpos
        kpoints
        parameters
        '''.split()
    assert(set(vxml.modeling.keys())==set(keys))

    parameters = obj(
        gga_compat      = True,
        i_constrained_m = 0,
        icorelevel      = 0,
        lberry          = False,
        ldau            = False,
        dos = obj(
            efermi          = 0.0,
            emax            = -10.0,
            emin            = 10.0,
            lorbit          = False,
            nedos           = 301,
            rwigs           = -1.0,
            ),
        electronic = obj(
            ediff           = 1e-06,
            enaug           = 358.966,
            enmax           = 230.283,
            eref            = 0.0,
            ialgo           = 38,
            irestart        = 0,
            iwavpr          = 11,
            nbands          = 96,
            nelect          = 130.0,
            nmin            = 0,
            nreboot         = 0,
            prec            = 'normal',
            turbo           = 0,
            electronic_convergence = obj(
                enini           = 230.283,
                nelm            = 60,
                nelmdl          = -5,
                nelmin          = 2,
                electronic_convergence_detail = obj(
                    deper           = 0.3,
                    ebreak          = 0.0,
                    ldiag           = True,
                    lsubrot         = True,
                    nrmm            = 4,
                    time            = 0.4,
                    weimin          = 0.001,
                    ),
                ),
            electronic_dipolcorrection = obj(
                dipol           = array([-100., -100., -100.],dtype=float),
                efield          = 0.0,
                epsilon         = 1.0,
                idipol          = 0,
                ldipol          = False,
                lmono           = False,
                ),
            electronic_exchange_correlation = obj(
                lasph           = False,
                lmetagga        = False,
                ),
            electronic_mixer = obj(
                amin            = 0.1,
                amix            = 0.4,
                amix_mag        = 1.6,
                bmix            = 1.0,
                bmix_mag        = 1.0,
                electronic_mixer_details = obj(
                    imix            = 4,
                    inimix          = 1,
                    maxmix          = -45,
                    mixpre          = 1,
                    mremove         = 5,
                    wc              = 100.0,
                    ),
                ),
            electronic_projectors = obj(
                lmaxmix         = 2,
                lmaxpaw         = -100,
                lreal           = False,
                nlspline        = False,
                ropt            = 0.0,
                ),
            electronic_smearing = obj(
                ismear          = 1,
                kgamma          = True,
                kspacing        = 0.5,
                sigma           = 0.2,
                ),
            electronic_spin = obj(
                ispin           = 1,
                lnoncollinear   = False,
                lsorbit         = False,
                lspiral         = False,
                lzeroz          = False,
                magmom          = array([1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.],dtype=float),
                nupdown         = -1.0,
                qspiral         = array([0.,0.,0.],dtype=float),
                saxis           = array([0.,0.,1.],dtype=float),
                ),
            electronic_startup = obj(
                icharg          = 2,
                iniwav          = 1,
                istart          = 0,
                ),
            ),
        electronic_exchange_correlation = obj(
            aexx            = 0.0,
            aggac           = 1.0,
            aggax           = 1.0,
            aldac           = 1.0,
            aldax           = 1.0,
            encut4o         = -1.0,
            evenonly        = False,
            exxoep          = 0,
            fourorbit       = 0,
            gga             = '--',
            hfalpha         = 0.0,
            hfscreen        = 0.0,
            hfscreenc       = 0.0,
            lfockaedft      = False,
            lhfcalc         = False,
            lhfone          = False,
            lmaxfock        = 0,
            lmodelhf        = False,
            lrhfcalc        = False,
            lsymgrad        = False,
            lthomas         = False,
            mcalpha         = 0.0,
            nbandsgwlow     = 0,
            nkredx          = 1,
            nkredy          = 1,
            nkredz          = 1,
            nmaxfockae      = 0,
            oddonly         = False,
            precfock        = '',
            shiftred        = False,
            voskown         = 0,
            ),
        general = obj(
            lcompat         = False,
            system          = array(['unknown','system'],dtype=str),
            ),
        grids = obj(
            addgrid         = False,
            ngx             = 24,
            ngxf            = 48,
            ngy             = 24,
            ngyf            = 48,
            ngz             = 64,
            ngzf            = 128,
            ),
        ionic = obj(
            ediffg          = -0.01,
            ibrion          = 1,
            isif            = 2,
            nfree           = 0,
            nsw             = 20,
            potim           = 0.5,
            pstress         = 0.0,
            scalee          = 1.0,
            smass           = -3.0,
            ),
        ionic_md = obj(
            apaco           = 16.0,
            kblock          = 20,
            nblock          = 1,
            npaco           = 256,
            tebeg           = 0.0001,
            teend           = 0.0001,
            ),
        linear_response_parameters = obj(
            cshift          = 0.1,
            deg_threshold   = 0.002,
            kinter          = 0,
            lepsilon        = False,
            lnabla          = False,
            lrpa            = False,
            lvel            = False,
            omegamax        = -1.0,
            ),
        miscellaneous = obj(
            darwinr         = 0.0,
            darwinv         = 1.0,
            idiot           = 3,
            lcorr           = True,
            lmusic          = False,
            pomass          = 195.08,
            ),
        model_gw = obj(
            model_alpha     = 1.0,
            model_eps0      = 24.58067119,
            model_gw        = 0,
            ),
        orbital_magnetization = obj(
            avecconst       = array([0.,0.,0.],dtype=float),
            lchimag         = False,
            lgauge          = True,
            lmagbloch       = False,
            lnicsall        = True,
            magatom         = 0,
            magdipol        = array([0.,0.,0.],dtype=float),
            magpos          = array([0.,0.,0.],dtype=float),
            nucind          = False,
            orbitalmag      = False,
            ),
        performance = obj(
            lasync          = False,
            lorbitalreal    = False,
            lplane          = True,
            lscaaware       = False,
            lscalapack      = False,
            lscalu          = False,
            nblk            = -1,
            npar            = 32,
            nsim            = 4,
            ),
        response_functions = obj(
            antires         = 0,
            cshift          = -0.1,
            dim             = 3,
            encutgw         = -2.0,
            encutgwsoft     = -2.0,
            encutlf         = -1.0,
            evenonlygw      = False,
            ibse            = 0,
            l2order         = False,
            ladder          = False,
            lfxc            = False,
            lfxceps         = False,
            lhartree        = True,
            lmaxmp2         = -1,
            lspectral       = False,
            ltcte           = False,
            ltete           = False,
            ltriplet        = False,
            lusew           = False,
            maxmem          = 1800,
            nbandsgw        = -1,
            nbandso         = -1,
            nbandsv         = -1,
            nelm            = 1,
            nkredlfx        = 1,
            nkredlfy        = 1,
            nkredlfz        = 1,
            nomega          = 0,
            nomegar         = 0,
            oddonlygw       = False,
            omegagrid       = 0,
            omegamax        = -30.0,
            omegatl         = -200.0,
            scissor         = 0.0,
            selfenergy      = False,
            telescope       = 0,
            ),
        symmetry = obj(
            isym            = 2,
            symprec         = 1e-05,
            ),
        vdw_dft = obj(
            luse_vdw        = False,
            param1          = 0.1234,
            param2          = 1.0,
            param3          = 0.0,
            zab_vdw         = -0.8491,
            ),
        writing = obj(
            lcharg          = True,
            lelf            = False,
            loptics         = False,
            lpard           = False,
            lvhar           = False,
            lvtot           = False,
            lwave           = True,
            nwrite          = 2,
            stm             = array([0.,0.,0.,0.,0.,0.,0.],dtype=float),
            ),
        )

    assert(object_eq(vxml.modeling.parameters,parameters))

#end def test_analyze
