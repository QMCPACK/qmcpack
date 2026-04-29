try:
    import pytest
    from . import NexusTestOrder
    pytestmark = pytest.mark.order(NexusTestOrder.PWSCF_INPUT)
    from ..generic import generic_settings
    generic_settings.raise_error = True
except ImportError:
    pass

from pathlib import Path
from ..testing import failed
from ..testing import divert_nexus_log,restore_nexus_log
from ..testing import object_eq,object_diff


associated_files = dict()

TEST_FILES = {
    "Cr_noncolin.in":         Path(__file__+"/../test_pwscf_input_files/Cr_noncolin.in").resolve(),
    "Fe_start_ns_eig.in":     Path(__file__+"/../test_pwscf_input_files/Fe_start_ns_eig.in").resolve(),
    "LiI_vc_relax.in":        Path(__file__+"/../test_pwscf_input_files/LiI_vc_relax.in").resolve(),
    "Ni_surface.in":          Path(__file__+"/../test_pwscf_input_files/Ni_surface.in").resolve(),
    "TiO2_band_structure.in": Path(__file__+"/../test_pwscf_input_files/TiO2_band_structure.in").resolve(),
    "TiO2_relax_freeze.in":   Path(__file__+"/../test_pwscf_input_files/TiO2_relax_freeze.in").resolve(),
    "VO2_M1_afm.in":          Path(__file__+"/../test_pwscf_input_files/VO2_M1_afm.in").resolve(),
    "VO2_M1_afm.xsf":         Path(__file__+"/../test_pwscf_input_files/VO2_M1_afm.xsf").resolve(),
    "WSe2_band_structure.in": Path(__file__+"/../test_pwscf_input_files/WSe2_band_structure.in").resolve(),
    "README": Path(__file__+"/../test_pwscf_input_files/README").resolve(),
}


def test_input(tmp_path):
    # imports
    import numpy as np
    from ..developer import obj
    from ..structure import read_structure
    from ..physical_system import generate_physical_system
    from ..pwscf_input import check_new_variables,check_section_classes
    from ..pwscf_input import PwscfInput,generate_pwscf_input

    tmp_dir = tmp_path / "test_pwscf_input_output"
    tmp_dir.mkdir()

    # divert logging function
    divert_nexus_log()

    # definitions
    def check_pw_same(pw1_,pw2_,l1='pw1',l2='pw2'):
        pw_same = object_eq(pw1_, pw2_, int_as_float=True, atol=5e-4)

        if not pw_same:
            d,d1,d2 = object_diff(pw1_,pw2_,full=True,int_as_float=True)
            diff = obj({l1:obj(d1),l2:obj(d2)})
            failed(str(diff))
        #end if
    #end def check_pw_same


    # test internal spec
    check_new_variables(exit=False)
    check_section_classes(exit=False)


    # test compose
    compositions = obj()

    # based on sample_inputs/Fe_start_ns_eig.in
    pw = PwscfInput()
    pw.control.set(
        calculation   = 'scf' ,
        restart_mode  = 'from_scratch' ,
        wf_collect    = True ,
        outdir        = './output' ,
        pseudo_dir    = '../pseudo/' ,
        prefix        = 'fe' ,
        etot_conv_thr = 1.0e-9 ,
        forc_conv_thr = 1.0e-6 ,
        tstress       = True ,
        tprnfor       = True ,
        )
    pw.system.set(
        ibrav           = 1,
        nat             = 2,
        ntyp            = 1,
        ecutwfc         = 100 ,
        ecutrho         = 300 ,
        nbnd            = 18,
        occupations     = 'smearing',
        degauss         = 0.0005 ,
        smearing        = 'methfessel-paxton' ,
        nspin           = 2 ,
        assume_isolated = 'martyna-tuckerman',
        lda_plus_u      = True ,
        )
    pw.system.set({
        'celldm(1)' : 15,
        'starting_magnetization(1)' : 0.9,
        'hubbard_u(1)' : 3.1,
        'starting_ns_eigenvalue(1,2,1)' : 0.0,
        'starting_ns_eigenvalue(2,2,1)' : 0.0476060,
        'starting_ns_eigenvalue(3,2,1)' : 0.0476060,
        'starting_ns_eigenvalue(4,2,1)' : 0.9654373,
        'starting_ns_eigenvalue(5,2,1)' : 0.9954307,
        })
    pw.electrons.set(
        conv_thr        = 1.0e-9 ,
        mixing_beta     = 0.7 ,
        diagonalization = 'david' ,
        mixing_fixed_ns = 500,
        )
    pw.atomic_species.set(
        atoms            = ['Fe'],
        masses           = obj(Fe=58.69000),
        pseudopotentials = obj(Fe='Fe.pbe-nd-rrkjus.UPF'),
        )
    pw.atomic_positions.set(
        specifier = 'angstrom',
        atoms     = ['Fe','Fe'],
        positions = np.array([
            [2.070000000,   0.000000000,   0.000000000],   
            [0.000000000,   0.000000000,   0.000000000], 
            ]),
        )
    pw.k_points.set(
        specifier = 'automatic',
        grid      = np.array((1,1,1)),
        shift     = np.array((1,1,1)),
        )

    compositions['Fe_start_ns_eig.in'] = pw


    # test read
    pwr = PwscfInput(TEST_FILES['Fe_start_ns_eig.in'])
    pwc = pw.copy()
    pwc.standardize_types()
    check_pw_same(pwc,pwr,'compose','read')


    # test write
    infile = tmp_dir / 'pwscf.in'
    pw.write(infile)
    pw2 = PwscfInput()
    pw2.read(infile)
    check_pw_same(pw2,pwr)


    # test read/write/read
    reads = obj()
    for input_file, file_path in TEST_FILES.items():
        if file_path.suffix != ".in":
            continue
        read_path = file_path
        write_path = tmp_dir / input_file
        pw = PwscfInput(read_path)
        pw.write(write_path)
        pw2 = PwscfInput(write_path)
        check_pw_same(pw,pw2,'read','write/read')
        reads[input_file] = pw
    #end for


    # test generate
    generations = obj()

    # based on sample_inputs/VO2_M1_afm.in
    infile      = 'VO2_M1_afm.in'
    struct_file = TEST_FILES['VO2_M1_afm.xsf']
    read_path   = TEST_FILES[infile]
    write_path  = tmp_dir / infile

    s = read_structure(struct_file)
    s.elem[0] = 'V1'
    s.elem[1] = 'V2'
    s.elem[2] = 'V1'
    s.elem[3] = 'V2'

    vo2 = generate_physical_system(
        structure = s,
        V1        = 13,
        V2        = 13,
        O         =  6,
        )

    pw = generate_pwscf_input(
        selector         = 'generic',
        calculation      = 'scf',
        disk_io          = 'low',
        verbosity        = 'high',
        wf_collect       = True,
        input_dft        = 'lda',
        hubbard_u        = obj(V1=3.5,V2=3.5),
        ecutwfc          = 350,
        bandfac          = 1.3,
        nosym            = True,
        occupations      = 'smearing',
        smearing         = 'fermi-dirac',
        degauss          = 0.0001,
        nspin            = 2,
        start_mag        = obj(V1=1.0,V2=-1.0),
        diagonalization  = 'david',
        conv_thr         = 1e-8,
        mixing_beta      = 0.2,
        electron_maxstep = 1000,
        system           = vo2,
        pseudos          = ['V.opt.upf','O.opt.upf'],
        kgrid            = (6,6,6),
        kshift           = (0,0,0),
        # added for reverse compatibility
        celldm           = {1:1.0},
        cell_option      = 'alat',
        positions_option = 'alat',
        )

    generations[infile] = pw

    pw.write(write_path)
    pw2 = PwscfInput(read_path)
    pw3 = PwscfInput(write_path)
    check_pw_same(pw2,pw3,'generate','read')

    # based on sample_inputs/Fe_start_ns_eig.in
    infile     = 'Fe_start_ns_eig.in'
    read_path  = TEST_FILES[infile]
    write_path = tmp_dir / infile

    pw = generate_pwscf_input(
        selector        = 'generic',
        calculation     = 'scf',
        restart_mode    = 'from_scratch',
        wf_collect      = True,
        outdir          = './output',
        pseudo_dir      = '../pseudo/',
        prefix          = 'fe',
        etot_conv_thr   = 1.0e-9,
        forc_conv_thr   = 1.0e-6,
        tstress         = True,
        tprnfor         = True,
        ibrav           = 1,
        nat             = 2,
        ntyp            = 1,
        ecutwfc         = 100,
        ecutrho         = 300,
        nbnd            = 18,
        occupations     = 'smearing',
        degauss         = 0.0005,
        smearing        = 'methfessel-paxton',
        nspin           = 2,
        assume_isolated = 'martyna-tuckerman',
        lda_plus_u      = True,
        conv_thr        = 1.0e-9,
        mixing_beta     = 0.7,
        diagonalization = 'david',
        mixing_fixed_ns = 500,
        mass            = obj(Fe=58.69000),
        pseudos         = ['Fe.pbe-nd-rrkjus.UPF'],
        elem            = ['Fe','Fe'],
        pos             = [[2.070000000, 0.000000000, 0.000000000],    
                           [0.000000000, 0.000000000, 0.000000000]],
        pos_specifier   = 'angstrom',
        kgrid           = np.array((1,1,1)),
        kshift          = np.array((1,1,1)),
        )
    pw.system.set({
        'celldm(1)' : 15,
        'starting_magnetization(1)' : 0.9,
        'hubbard_u(1)' : 3.1,
        'starting_ns_eigenvalue(1,2,1)' : 0.0,
        'starting_ns_eigenvalue(2,2,1)' : 0.0476060,
        'starting_ns_eigenvalue(3,2,1)' : 0.0476060,
        'starting_ns_eigenvalue(4,2,1)' : 0.9654373,
        'starting_ns_eigenvalue(5,2,1)' : 0.9954307,
        })

    generations[infile] = pw

    pw2 = compositions[infile]
    check_pw_same(pw,pw2,'generate','compose')
    pw.write(write_path)
    pw3 = PwscfInput(write_path)
    pw4 = reads[infile]
    check_pw_same(pw3,pw4,'generate','read')


    # based on sample_inputs/Fe_start_ns_eig.in
    #  variant that uses direct pwscf array input
    pw = generate_pwscf_input(
        selector        = 'generic',
        calculation     = 'scf',
        restart_mode    = 'from_scratch',
        wf_collect      = True,
        outdir          = './output',
        pseudo_dir      = '../pseudo/',
        prefix          = 'fe',
        etot_conv_thr   = 1.0e-9,
        forc_conv_thr   = 1.0e-6,
        tstress         = True,
        tprnfor         = True,
        ibrav           = 1,
        nat             = 2,
        ntyp            = 1,
        ecutwfc         = 100,
        ecutrho         = 300,
        nbnd            = 18,
        occupations     = 'smearing',
        degauss         = 0.0005,
        smearing        = 'methfessel-paxton',
        nspin           = 2,
        assume_isolated = 'martyna-tuckerman',
        lda_plus_u      = True,
        conv_thr        = 1.0e-9,
        mixing_beta     = 0.7,
        diagonalization = 'david',
        mixing_fixed_ns = 500,
        mass            = obj(Fe=58.69000),
        pseudos         = ['Fe.pbe-nd-rrkjus.UPF'],
        elem            = ['Fe','Fe'],
        pos             = [[2.070000000, 0.000000000, 0.000000000],    
                           [0.000000000, 0.000000000, 0.000000000]],
        pos_specifier   = 'angstrom',
        kgrid           = np.array((1,1,1)),
        kshift          = np.array((1,1,1)),
        starting_ns_eigenvalue = {(1,2,1) : 0.0,
                                  (2,2,1) : 0.0476060,
                                  (3,2,1) : 0.0476060,
                                  (4,2,1) : 0.9654373,
                                  (5,2,1) : 0.9954307,},
        celldm                 = {1 : 15 },
        starting_magnetization = {1 : 0.9},
        hubbard_u              = {1 : 3.1},
        )

    pwg = pw.copy()
    pwg.standardize_types()

    generations[infile] = pw

    pw2 = compositions[infile].copy()
    pw2.standardize_types()
    check_pw_same(pwg,pw2,'generate','compose')
    pw3 = reads[infile]
    check_pw_same(pwg,pw3,'generate','read')
    pw.write(write_path)
    pw4 = PwscfInput(write_path)
    check_pw_same(pwg,pw3,'generate','write')


    # restore logging function
    restore_nexus_log()
#end def test_input
