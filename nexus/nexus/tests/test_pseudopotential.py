import pytest
from copy import deepcopy
from . import NexusTestOrder
pytestmark = pytest.mark.order(NexusTestOrder.PSEUDOPOTENTIAL)


from . import TEST_DIR
from ..testing import value_eq, object_eq


TEST_FILES = {
    "C.BFD.gms": TEST_DIR / "test_pseudopotential_files/C.BFD.gms",
    "C.BFD.upf": TEST_DIR / "../examples/qmcpack/pseudopotentials/C.BFD.upf",
    "C.BFD.xml": TEST_DIR / "../examples/qmcpack/pseudopotentials/C.BFD.xml",
    }

for file in TEST_FILES.values():
    assert(file.exists()), f"Test file not found! {file}"


def test_pseudopotential_classes(tmp_path):
    import numpy as np
    from ..developer import to_obj
    from ..pseudopotential import SemilocalPP
    from ..pseudopotential import GaussianPP
    from ..pseudopotential import QmcpackPP
    from ..pseudopotential import CasinoPP

    # empty initialization
    SemilocalPP()
    GaussianPP()
    QmcpackPP()
    CasinoPP()

    qpp = QmcpackPP(TEST_FILES['C.BFD.xml'])

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

    qpp_fake = deepcopy(qpp)
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
           4.00000000000000e+00   4.00000000000000e+00   4.00000000000000e+00
           4.00000000000000e+00   4.00000000000000e+00   4.00000000000000e+00
        </data>
      </radfunc>
    </vps>
    <vps principal-n="0" l="p" spin="-1" cutoff="10.0" occupation="unknown">
      <radfunc>
        <grid type="linear" units="bohr" ri="0.0" rf="10.0" npts="6"/>
        <data>
           4.00000000000000e+00   4.00000000000000e+00   4.00000000000000e+00
           4.00000000000000e+00   4.00000000000000e+00   4.00000000000000e+00
        </data>
      </radfunc>
    </vps>
  </semilocal>
</pseudo>'''

    qtext = qpp_fake.write_qmcpack()
    assert(qtext.strip()==qtext_ref.strip())

    # Read legacy QMCpack files in which adjacent L2 values ran together
    # because a value filled the allotted output field width.
    l2_text = '''  <L2 units="hartree" format="r*V" cutoff="10.0">
    <radfunc>
      <grid type="linear" units="bohr" ri="0.0" rf="10.0" npts="3"/>
      <data>
         1.00000000000000e+00-2.00000000000000e+00 3.00000000000000e+00
      </data>
    </radfunc>
  </L2>
'''
    qtext_legacy_l2 = qtext.replace('  <semilocal',l2_text+'  <semilocal')
    legacy_l2_file = tmp_path / 'legacy_l2.qmcpack'
    legacy_l2_file.write_text(qtext_legacy_l2)
    qpp_legacy_l2 = QmcpackPP(legacy_l2_file)
    assert(qpp_legacy_l2.has_L2())
    assert(value_eq(qpp_legacy_l2.components.L2,np.array([1.,-2.,3.])))

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
    gpp = GaussianPP(TEST_FILES['C.BFD.gms'],format='gamess')
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
    gamess_file = tmp_path / 'C.BFD.gamess'
    gpp.write(gamess_file,format='gamess')

    gaussian_file = tmp_path / 'C.BFD.gaussian'
    gpp.write(gaussian_file,format='gaussian')

    qmcpack_file = tmp_path / 'C.BFD.qmcpack'
    gpp.write(qmcpack_file,format='qmcpack')

    casino_file = tmp_path / 'C.BFD.casino'
    gpp.write(casino_file,format='casino')


    gpp_gamess = GaussianPP(gamess_file,format='gamess')
    assert(object_eq(gpp_gamess,gpp))

    gpp_gaussian = GaussianPP(gaussian_file,format='gaussian')
    assert(object_eq(gpp_gaussian,gpp))

    qpp_qmcpack = QmcpackPP(qmcpack_file)
    assert(object_eq(qpp_qmcpack,qpp,int_as_float=True,atol=1e-12))


    # tests for CasinoPP
    cpp = CasinoPP(casino_file)

    qo = to_obj(qpp)
    co = to_obj(cpp)
    del qo.rmin
    del qo.rmax
    assert(object_eq(co,qo,atol=1e-12))

    qmcpack_from_casino_file = tmp_path / 'C.BFD.qmcpack_from_casino'
    cpp.write_qmcpack(qmcpack_from_casino_file)

    qpp_casino = QmcpackPP(qmcpack_from_casino_file)
    assert(object_eq(qpp_casino,qpp,int_as_float=True,atol=1e-12))

#end def test_pseudopotential_classes
