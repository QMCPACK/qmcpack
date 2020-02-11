
import testing
from testing import value_eq,object_eq
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



def test_ppset():
    from generic import obj
    from pseudopotential import ppset

    ppset_ref = obj(
        pseudos = obj(
            bfd = obj(
                gamess  = obj(C='C.BFD.gms'),
                pwscf   = obj(C='C.BFD.upf'),
                qmcpack = obj(C='C.BFD.xml'),
                ),
            ),
        )

    ppset(
        label   = 'bfd',
        gamess  = ['C.BFD.gms'],
        pwscf   = ['C.BFD.upf'],
        qmcpack = ['C.BFD.xml'],
        )

    o = ppset.to_obj()
    assert(object_eq(o,ppset_ref))

    assert(ppset.supports_code('pwscf'))
    assert(ppset.supports_code('gamess'))
    assert(ppset.supports_code('vasp'))
    assert(ppset.supports_code('qmcpack'))

    assert(ppset.has_set('bfd'))


    # need to add test for get() method
    #   depends on PhysicalSystem
#end def test_ppset



def test_pseudopotential_classes():
    import os
    import numpy as np
    from pseudopotential import SemilocalPP
    from pseudopotential import GaussianPP
    from pseudopotential import QmcpackPP
    from pseudopotential import CasinoPP

    tpath = testing.setup_unit_test_output_directory('pseudopotential','test_pseudopotential_classes')

    files = get_files()

    # empty initialization
    SemilocalPP()
    GaussianPP()
    QmcpackPP()
    CasinoPP()


    qpp = QmcpackPP(files['C.BFD.xml'])

    # SemilocalPP attributes/methods
    assert(qpp.name is None)
    assert(qpp.rcut is None)
    assert(qpp.lmax==1)
    assert(qpp.local=='p')

    assert(qpp.has_component('s'))
    assert(qpp.has_component('p'))

    assert(isinstance(qpp.get_component('s'),np.ndarray))
    assert(isinstance(qpp.get_component('p'),np.ndarray))

    assert(qpp.has_local())
    assert(qpp.has_nonlocal())
    assert(not qpp.has_L2())

    vnl = qpp.get_nonlocal()

    rc = qpp.find_rcut()
    assert(value_eq(rc,1.705,atol=1e-3))

    # below follows by virtue of being numeric
    qpp.assert_numeric('some location')
    
    vcomp = qpp.components

    vloc = qpp.evaluate_local(rpow=1)
    assert(value_eq(vloc,vcomp.p))

    vnonloc = qpp.evaluate_nonlocal(l='s',rpow=1)
    assert(value_eq(vnonloc,vcomp.s))

    vs = qpp.evaluate_channel(l='s',rpow=1)
    assert(value_eq(vs,vcomp.s+vcomp.p))

    vp = qpp.evaluate_channel(l='p',rpow=1)
    assert(value_eq(vp,vcomp.p))

    r,vsn = qpp.numeric_channel(l='s',rpow=1)
    r,vpn = qpp.numeric_channel(l='p',rpow=1)
    assert(value_eq(r,qpp.r))
    assert(value_eq(vsn,vs))
    assert(value_eq(vpn,vp))

    # QmcpackPP attributes/methods
    assert(qpp.numeric)
    assert(qpp.Zcore==2)
    assert(qpp.Zval==4)
    assert(qpp.core=='He')
    assert(qpp.element=='C')
    assert(value_eq(float(qpp.rmin),0.))
    assert(value_eq(qpp.rmax,10.))
    assert(value_eq(qpp.r.min(),0.))
    assert(value_eq(qpp.r.max(),10.))

    assert(value_eq(qpp.v_at_zero('s'),22.551641791033372))
    assert(value_eq(qpp.v_at_zero('p'),-19.175372435022126))

    qpp_fake = qpp.copy()
    r = np.linspace(0,10,6)
    vloc = 0*r + qpp.Zval
    vnl  = 0*r
    qpp_fake.r = r
    qpp_fake.components.s = vnl
    qpp_fake.components.p = vloc

    qtext_ref = '''<?xml version="1.0" encoding="UTF-8"?>
<pseudo version="0.5">
  <header symbol="C" atomic-number="6" zval="4" relativistic="unknown" 
   polarized="unknown" creator="Nexus" flavor="unknown" 
   core-corrections="unknown" xc-functional-type="unknown" 
   xc-functional-parametrization="unknown"/>
  <grid type="linear" units="bohr" ri="0.0" rf="10.0" npts="6"/>
  <semilocal units="hartree" format="r*V" npots-down="2" npots-up="0" l-local="1">
    <vps principal-n="0" l="s" spin="-1" cutoff="10.0" occupation="unknown">
      <radfunc>
        <grid type="linear" units="bohr" ri="0.0" rf="10.0" npts="6"/>
        <data>
          4.00000000000000e+00  4.00000000000000e+00  4.00000000000000e+00
          4.00000000000000e+00  4.00000000000000e+00  4.00000000000000e+00
        </data>
      </radfunc>
    </vps>
    <vps principal-n="0" l="p" spin="-1" cutoff="10.0" occupation="unknown">
      <radfunc>
        <grid type="linear" units="bohr" ri="0.0" rf="10.0" npts="6"/>
        <data>
          4.00000000000000e+00  4.00000000000000e+00  4.00000000000000e+00
          4.00000000000000e+00  4.00000000000000e+00  4.00000000000000e+00
        </data>
      </radfunc>
    </vps>
  </semilocal>
</pseudo>'''

    qtext = qpp_fake.write_qmcpack()
    assert(qtext.strip()==qtext_ref.strip())

    ctext_ref = '''C pseudopotential converted by Nexus
Atomic number and pseudo-charge
  6 4.0
Energy units (rydberg/hartree/ev):
  hartree
Angular momentum of local component (0=s,1=p,2=d..)
  1
NLRULE override (1) VMC/DMC (2) config gen (0 ==> input/default value)
  0 0
Number of grid points
  6
R(i) in atomic units
  0.00000000000000e+00
  2.00000000000000e+00
  4.00000000000000e+00
  6.00000000000000e+00
  8.00000000000000e+00
  1.00000000000000e+01
r*potential (L=0) in Ha
  4.00000000000000e+00
  4.00000000000000e+00
  4.00000000000000e+00
  4.00000000000000e+00
  4.00000000000000e+00
  4.00000000000000e+00
r*potential (L=1) in Ha
  4.00000000000000e+00
  4.00000000000000e+00
  4.00000000000000e+00
  4.00000000000000e+00
  4.00000000000000e+00
  4.00000000000000e+00'''

    ctext = qpp_fake.write_casino()
    assert(ctext.strip()==ctext_ref.strip())

    
    # tests for GaussianPP
    gpp = GaussianPP(files['C.BFD.gms'],format='gamess')
    assert(gpp.Zcore   == 2   )
    assert(gpp.Zval    == 4   )
    assert(gpp.core    == 'He')
    assert(gpp.element == 'C' )
    assert(gpp.lmax    == 1   )
    assert(gpp.local   == 'p' )
    assert(gpp.name is None)
    assert(value_eq(gpp.rcut,1.7053,atol=1e-3))

    assert(len(gpp.basis)==20)

    nterms_ref = [9,1,1,1,1,9,1,1,1,1,1,1,1,1,1,1,1,1,1,1]
    nterms = []
    for n in range(len(gpp.basis)):
        nterms.append(len(gpp.basis[n].terms))
    #end for
    assert(nterms==nterms_ref)

    assert(value_eq(gpp.basis[5].terms[4].coeff,0.289868))

    assert(len(gpp.components.s)==1)
    assert(len(gpp.components.p)==3)

    assert(value_eq(gpp.components.p[1].expon,4.48361888))

    # check cross-format write/read
    gamess_file = os.path.join(tpath,'C.BFD.gamess')
    gpp.write(gamess_file,format='gamess')

    gaussian_file = os.path.join(tpath,'C.BFD.gaussian')
    gpp.write(gaussian_file,format='gaussian')

    qmcpack_file = os.path.join(tpath,'C.BFD.qmcpack')
    gpp.write(qmcpack_file,format='qmcpack')

    casino_file = os.path.join(tpath,'C.BFD.casino')
    gpp.write(casino_file,format='casino')


    gpp_gamess = GaussianPP(gamess_file,format='gamess')
    assert(object_eq(gpp_gamess,gpp))

    gpp_gaussian = GaussianPP(gaussian_file,format='gaussian')
    assert(object_eq(gpp_gaussian,gpp))

    qpp_qmcpack = QmcpackPP(qmcpack_file)
    assert(object_eq(qpp_qmcpack,qpp,int_as_float=True,atol=1e-12))


    # tests for CasinoPP
    cpp = CasinoPP(casino_file)

    qo = qpp.to_obj()
    co = cpp.to_obj()
    del qo.rmin
    del qo.rmax
    assert(object_eq(co,qo,atol=1e-12))

    qmcpack_from_casino_file = os.path.join(tpath,'C.BFD.qmcpack_from_casino')
    cpp.write_qmcpack(qmcpack_from_casino_file)

    qpp_casino = QmcpackPP(qmcpack_from_casino_file)
    assert(object_eq(qpp_casino,qpp,int_as_float=True,atol=1e-12))

#end def test_pseudopotential_classes

