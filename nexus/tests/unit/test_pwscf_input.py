
import testing
from testing import failed
from testing import divert_nexus_log,restore_nexus_log
from testing import value_eq,object_eq


associated_files = dict()

input_files = [
    'Cr_noncolin.in',
    'Fe_start_ns_eig.in',
    'LiI_vc_relax.in',
    'Ni_surface.in',
    'TiO2_band_structure.in',
    'TiO2_relax_freeze.in',
    'VO2_M1_afm.in',
    'WSe2_band_structure.in',
    ]

structure_files = [
    'VO2_M1_afm.xsf',
    ]

other_files = [
    'README',
    ]


def get_files():
    return testing.collect_unit_test_file_paths('pwscf_input',associated_files)
#end def get_files



def test_files():
    filenames = input_files + structure_files + other_files
    files = get_files()
    assert(set(files.keys())==set(filenames))
#end def test_files



def test_import():
    from pwscf_input import check_new_variables,check_section_classes
    from pwscf_input import PwscfInput,generate_pwscf_input
#end def test_import



def test_input():
    # imports
    import os
    import numpy as np
    import pwscf_input as pwi
    from generic import obj
    from structure import read_structure
    from physical_system import generate_physical_system
    from pwscf_input import check_new_variables,check_section_classes
    from pwscf_input import PwscfInput,generate_pwscf_input
 
    # directories
    tpath = testing.setup_unit_test_output_directory('pwscf_input','test_input')

    # files
    files = get_files()

    # divert logging function
    divert_nexus_log()

    # definitions
    def check_pw_same(pw1_,pw2_,l1='pw1',l2='pw2'):
        pw_same = object_eq(pw1_,pw2_,int_as_float=True)
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
    pwr = PwscfInput(files['Fe_start_ns_eig.in'])
    pwc = pw.copy()
    pwc.standardize_types()
    check_pw_same(pwc,pwr,'compose','read')


    # test write
    infile = os.path.join(tpath,'pwscf.in')
    pw.write(infile)
    pw2 = PwscfInput()
    pw2.read(infile)
    check_pw_same(pw2,pwr)


    # test read/write/read
    reads = obj()
    for infile in input_files:
        read_path = files[infile]
        write_path = os.path.join(tpath,infile)
        if os.path.exists(write_path):
            os.remove(write_path)
        #end if
        pw = PwscfInput(read_path)
        pw.write(write_path)
        pw2 = PwscfInput(write_path)
        check_pw_same(pw,pw2,'read','write/read')
        reads[infile] = pw
    #end for


    # test generate
    generations = obj()

    # based on sample_inputs/VO2_M1_afm.in
    infile      = 'VO2_M1_afm.in'
    struct_file = files['VO2_M1_afm.xsf']
    read_path   = files[infile]
    write_path  = os.path.join(tpath,infile)

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

    if os.path.exists(write_path):
        os.remove(write_path)
    #end if
    pw.write(write_path)
    pw2 = PwscfInput(read_path)
    pw3 = PwscfInput(write_path)
    check_pw_same(pw2,pw3,'generate','read')

    # based on sample_inputs/Fe_start_ns_eig.in
    infile     = 'Fe_start_ns_eig.in'
    read_path  = files[infile]
    write_path = os.path.join(tpath,infile)

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
    if os.path.exists(write_path):
        os.remove(write_path)
    #end if
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
    if os.path.exists(write_path):
        os.remove(write_path)
    #end if
    pw.write(write_path)
    pw4 = PwscfInput(write_path)
    check_pw_same(pwg,pw3,'generate','write')


    # restore logging function
    restore_nexus_log()
#end def test_input
