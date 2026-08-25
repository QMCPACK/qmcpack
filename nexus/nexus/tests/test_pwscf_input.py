import pytest
from copy import deepcopy
from . import NexusTestOrder
pytestmark = pytest.mark.order(NexusTestOrder.PWSCF_INPUT)


from . import isolate_nexus_core, register_pseudo_files, TEST_DIR
from ..testing import failed
from ..testing import object_eq,object_diff


TEST_FILES = {
    "Cr_noncolin.in":         TEST_DIR / "test_pwscf_input_files/Cr_noncolin.in",
    "Fe_start_ns_eig.in":     TEST_DIR / "test_pwscf_input_files/Fe_start_ns_eig.in",
    "LiI_vc_relax.in":        TEST_DIR / "test_pwscf_input_files/LiI_vc_relax.in",
    "Ni_surface.in":          TEST_DIR / "test_pwscf_input_files/Ni_surface.in",
    "nexus_argon_scf.in":     TEST_DIR / "test_pwscf_input_files/nexus_argon_scf.in",
    "nexus_h2_relax.in":      TEST_DIR / "test_pwscf_input_files/nexus_h2_relax.in",
    "nexus_argon_bands.in":  TEST_DIR / "test_pwscf_input_files/nexus_argon_bands.in",
    "TiO2_band_structure.in": TEST_DIR / "test_pwscf_input_files/TiO2_band_structure.in",
    "TiO2_relax_freeze.in":   TEST_DIR / "test_pwscf_input_files/TiO2_relax_freeze.in",
    "VO2_M1_afm.in":          TEST_DIR / "test_pwscf_input_files/VO2_M1_afm.in",
    "VO2_M1_afm.xsf":         TEST_DIR / "test_pwscf_input_files/VO2_M1_afm.xsf",
    "WSe2_band_structure.in": TEST_DIR / "test_pwscf_input_files/WSe2_band_structure.in",
    "README": TEST_DIR / "test_pwscf_input_files/README",
    }

for file in TEST_FILES.values():
    assert(file.exists()), f"Test file not found! {file}"

@isolate_nexus_core
def test_input(tmp_path):
    register_pseudo_files([
        'V.opt.upf','O.opt.upf','Fe.pbe-nd-rrkjus.UPF'
        ])
    # imports
    import numpy as np
    from ..developer import obj
    from ..structure import read_structure
    from ..physical_system import generate_physical_system
    from ..pwscf_input import check_section_classes
    from ..pwscf_input import PwscfInput,generate_pwscf_input

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
    check_section_classes(exit=False)


    # test compose
    compositions = obj()

    # based on sample_inputs/Fe_start_ns_eig.in
    pw = PwscfInput()
    pw.control.update(
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
    pw.system.update(
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
    pw.system.update({
        'celldm(1)' : 15,
        'starting_magnetization(1)' : 0.9,
        'hubbard_u(1)' : 3.1,
        'starting_ns_eigenvalue(1,2,1)' : 0.0,
        'starting_ns_eigenvalue(2,2,1)' : 0.0476060,
        'starting_ns_eigenvalue(3,2,1)' : 0.0476060,
        'starting_ns_eigenvalue(4,2,1)' : 0.9654373,
        'starting_ns_eigenvalue(5,2,1)' : 0.9954307,
        })
    pw.electrons.update(
        conv_thr        = 1.0e-9 ,
        mixing_beta     = 0.7 ,
        diagonalization = 'david' ,
        mixing_fixed_ns = 500,
        )
    pw.atomic_species.update(
        atoms            = ['Fe'],
        masses           = obj(Fe=58.69000),
        pseudopotentials = obj(Fe='Fe.pbe-nd-rrkjus.UPF'),
        )
    pw.atomic_positions.update(
        specifier = 'angstrom',
        atoms     = ['Fe','Fe'],
        positions = np.array([
            [2.070000000,   0.000000000,   0.000000000],   
            [0.000000000,   0.000000000,   0.000000000], 
            ]),
        )
    pw.k_points.update(
        specifier = 'automatic',
        grid      = np.array((1,1,1)),
        shift     = np.array((1,1,1)),
        )

    compositions['Fe_start_ns_eig.in'] = pw

    # Nexus-authored representative scf input
    pw = PwscfInput()
    pw.control.update(
        calculation = 'scf',
        prefix       = 'nexus_argon',
        outdir       = './nexus_test_tmp',
        pseudo_dir   = './pseudo',
        )
    pw.system.update(
        ibrav       = 1,
        nat         = 1,
        ntyp        = 1,
        nspin       = 1,
        ecutwfc     = 32.0,
        occupations = 'fixed',
        )
    pw.system['celldm(1)'] = 9.25
    pw.electrons.update(
        conv_thr        = 2.0e-9,
        diagonalization = 'david',
        )
    pw.atomic_species.update(
        atoms            = ['Ar'],
        masses           = obj(Ar=39.948),
        pseudopotentials = obj(Ar='Ar.nexus-test.UPF'),
        )
    pw.atomic_positions.update(
        specifier = 'crystal',
        atoms     = ['Ar'],
        positions = np.array([[0.0,0.0,0.0]]),
        )
    pw.k_points.update(
        specifier = 'automatic',
        grid      = np.array((5,5,5)),
        shift     = np.array((1,1,1)),
        )
    compositions['nexus_argon_scf.in'] = pw

    # Nexus-authored representative relax input
    pw = PwscfInput('ions')
    pw.control.update(
        calculation   = 'relax',
        prefix        = 'nexus_h2',
        outdir        = './nexus_test_tmp',
        pseudo_dir    = './pseudo',
        forc_conv_thr = 3.0e-5,
        nstep         = 40,
        )
    pw.system.update(
        ibrav           = 1,
        nat             = 2,
        ntyp            = 1,
        nspin           = 1,
        ecutwfc         = 28.0,
        assume_isolated = 'martyna-tuckerman',
        )
    pw.system['celldm(1)'] = 18.0
    pw.electrons.update(
        conv_thr    = 4.0e-10,
        mixing_beta = 0.55,
        )
    pw.ions.ion_dynamics = 'bfgs'
    pw.atomic_species.update(
        atoms            = ['H'],
        masses           = obj(H=1.00794),
        pseudopotentials = obj(H='H.nexus-test.UPF'),
        )
    pw.atomic_positions.update(
        specifier        = 'angstrom',
        atoms            = ['H','H'],
        positions        = np.array([
            [4.50,4.50,4.09],
            [4.50,4.50,4.91],
            ]),
        relax_directions = np.array([
            [0,0,1],
            [0,0,1],
            ]),
        )
    pw.k_points.specifier = 'gamma'
    compositions['nexus_h2_relax.in'] = pw

    # Nexus-authored representative bands input
    pw = PwscfInput()
    pw.control.update(
        calculation = 'bands',
        prefix       = 'nexus_argon',
        outdir       = './nexus_test_tmp',
        pseudo_dir   = './pseudo',
        )
    pw.system.update(
        ibrav       = 1,
        nat         = 1,
        ntyp        = 1,
        nspin       = 1,
        ecutwfc     = 32.0,
        nbnd        = 8,
        occupations = 'fixed',
        )
    pw.system['celldm(1)'] = 9.25
    pw.electrons.update(
        conv_thr        = 2.0e-9,
        diagonalization = 'david',
        )
    pw.atomic_species.update(
        atoms            = ['Ar'],
        masses           = obj(Ar=39.948),
        pseudopotentials = obj(Ar='Ar.nexus-test.UPF'),
        )
    pw.atomic_positions.update(
        specifier = 'crystal',
        atoms     = ['Ar'],
        positions = np.array([[0.0,0.0,0.0]]),
        )
    pw.k_points.update(
        specifier = 'crystal_b',
        nkpoints  = 4,
        kpoints   = np.array([
            [0.0,0.0,0.0],
            [0.5,0.0,0.0],
            [0.5,0.5,0.0],
            [0.0,0.0,0.0],
            ]),
        weights   = np.array((16,16,16,1)),
        )
    compositions['nexus_argon_bands.in'] = pw


    # test read
    for input_file,pw in compositions.items():
        pwr = PwscfInput(TEST_FILES[input_file])
        pwc = deepcopy(pw)
        pwc.standardize_types()
        check_pw_same(pwc,pwr,'compose','read')
    #end for


    # test write
    for input_file,pw in compositions.items():
        write_path = tmp_path / ('composed_'+input_file)
        pw.write(write_path)
        pw2 = PwscfInput()
        pw2.read(write_path)
        pwr = PwscfInput(TEST_FILES[input_file])
        check_pw_same(pw2,pwr,'write/read','read')
    #end for


    # test read/write/read
    reads = obj()
    for input_file, file_path in TEST_FILES.items():
        if file_path.suffix != ".in":
            continue
        read_path = file_path
        write_path = tmp_path / input_file
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
    write_path  = tmp_path / infile

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
    write_path = tmp_path / infile

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
    pw.system.update({
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

    pwg = deepcopy(pw)
    pwg.standardize_types()

    generations[infile] = pw

    pw2 = deepcopy(compositions[infile])
    pw2.standardize_types()
    check_pw_same(pwg,pw2,'generate','compose')
    pw3 = reads[infile]
    check_pw_same(pwg,pw3,'generate','read')
    pw.write(write_path)
    pw4 = PwscfInput(write_path)
    check_pw_same(pwg,pw3,'generate','write')

    # Nexus-authored representative scf input
    infile     = 'nexus_argon_scf.in'
    write_path = tmp_path / infile
    pw = generate_pwscf_input(
        selector        = 'generic',
        calculation     = 'scf',
        prefix          = 'nexus_argon',
        outdir          = './nexus_test_tmp',
        pseudo_dir      = './pseudo',
        ibrav           = 1,
        celldm          = {1:9.25},
        nat             = 1,
        ntyp            = 1,
        nspin           = 1,
        ecutwfc         = 32.0,
        occupations     = 'fixed',
        conv_thr        = 2.0e-9,
        diagonalization = 'david',
        mass            = obj(Ar=39.948),
        pseudos         = ['Ar.nexus-test.UPF'],
        elem            = ['Ar'],
        pos             = [[0.0,0.0,0.0]],
        pos_specifier   = 'crystal',
        kgrid           = (5,5,5),
        kshift          = (1,1,1),
        )
    generations[infile] = pw
    pw.write(write_path)
    pw2 = PwscfInput(write_path)
    check_pw_same(pw2,reads[infile],'generate','read')

    # Nexus-authored representative relax input
    infile     = 'nexus_h2_relax.in'
    write_path = tmp_path / infile
    pw = generate_pwscf_input(
        selector        = 'generic',
        calculation     = 'relax',
        prefix          = 'nexus_h2',
        outdir          = './nexus_test_tmp',
        pseudo_dir      = './pseudo',
        forc_conv_thr   = 3.0e-5,
        nstep           = 40,
        ibrav           = 1,
        celldm          = {1:18.0},
        nat             = 2,
        ntyp            = 1,
        nspin           = 1,
        ecutwfc         = 28.0,
        assume_isolated = 'martyna-tuckerman',
        conv_thr        = 4.0e-10,
        mixing_beta     = 0.55,
        ion_dynamics    = 'bfgs',
        mass            = obj(H=1.00794),
        pseudos         = ['H.nexus-test.UPF'],
        elem            = ['H','H'],
        pos             = [[4.50,4.50,4.09],
                           [4.50,4.50,4.91]],
        pos_specifier   = 'angstrom',
        kgrid           = (1,1,1),
        kshift          = (0,0,0),
        )
    pw.atomic_positions.relax_directions = np.array([
        [0,0,1],
        [0,0,1],
        ])
    generations[infile] = pw
    pw.write(write_path)
    pw2 = PwscfInput(write_path)
    check_pw_same(pw2,reads[infile],'generate','read')

    # Nexus-authored representative bands input.  The explicit path is set on
    # the generated input because the generic interface only accepts regular
    # k-point grids directly.
    infile     = 'nexus_argon_bands.in'
    write_path = tmp_path / infile
    pw = generate_pwscf_input(
        selector        = 'generic',
        calculation     = 'bands',
        prefix          = 'nexus_argon',
        outdir          = './nexus_test_tmp',
        pseudo_dir      = './pseudo',
        ibrav           = 1,
        celldm          = {1:9.25},
        nat             = 1,
        ntyp            = 1,
        nspin           = 1,
        ecutwfc         = 32.0,
        nbnd            = 8,
        occupations     = 'fixed',
        conv_thr        = 2.0e-9,
        diagonalization = 'david',
        mass            = obj(Ar=39.948),
        pseudos         = ['Ar.nexus-test.UPF'],
        elem            = ['Ar'],
        pos             = [[0.0,0.0,0.0]],
        pos_specifier   = 'crystal',
        kgrid           = (1,1,1),
        kshift          = (0,0,0),
        )
    pw.k_points.update(
        specifier = 'crystal_b',
        nkpoints  = 4,
        kpoints   = np.array([
            [0.0,0.0,0.0],
            [0.5,0.0,0.0],
            [0.5,0.5,0.0],
            [0.0,0.0,0.0],
            ]),
        weights   = np.array((16,16,16,1)),
        )
    generations[infile] = pw
    pw.write(write_path)
    pw2 = PwscfInput(write_path)
    check_pw_same(pw2,reads[infile],'generate','read')
#end def test_input
