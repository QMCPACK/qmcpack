import pytest
from . import NexusTestOrder
pytestmark = pytest.mark.order(NexusTestOrder.PSEUDOPOTENTIAL)

from ..generic import generic_settings, NexusUserWarning
generic_settings.raise_error = True

from pathlib import Path
from typing import Literal
from . import isolate_nexus_core, TEST_DIR
from ..testing import value_eq,object_eq
from nexus.pseudopotential import read_upf_z_valence, read_xml_z_valence, read_potcar_z_valence
from nexus.pseudopotential import PseudoSet


TEST_FILES = {
    "C.BFD.gms": TEST_DIR / "test_pseudopotential_files/C.BFD.gms",
    "C.BFD.upf": TEST_DIR / "../examples/qmcpack/pseudopotentials/C.BFD.upf",
    "C.BFD.xml": TEST_DIR / "../examples/qmcpack/pseudopotentials/C.BFD.xml",
    }

for file in TEST_FILES.values():
    assert(file.exists()), f"Test file not found! {file}"


def setup_psps(
    test_dir: Path,
    code: Literal["espresso", "gamess", "vasp", "qmcpack"],
    ) -> tuple[Path, list[Path], dict[str, Path]]:
    """Take a test's temp directory and populate with dummy pseudopotential files."""
    match code:
        case "espresso":
            file_ext = "upf"
        case "qmcpack":
            file_ext = "xml"
        case "gamess":
            file_ext = "gms"
        case "vasp":
            file_ext = "potcar"
        case _:
            raise pytest.UsageError(
                "Invalid call to `setup_for_pseudoset()`!\n"
                f"Code supplied is {code}, but must be one of: {', '.join(PseudoSet.known_codes)}"
                )

    pseudo_names = (
        f"C.BFD.{file_ext}",
        f"H.BFD.{file_ext}",
        f"O.BFD.{file_ext}",
        )

    psp_dir = test_dir / f"{file_ext}_pseudos"
    psp_dir.mkdir()
    assert psp_dir.exists(), "Failed to create pseudo directory!"

    pseudo_list = []
    for psp in pseudo_names:
        pseudo = (psp_dir / psp).resolve()
        pseudo.touch()
        assert pseudo.exists(), "Failed to create pseudo file!"
        pseudo_list.append(pseudo)

    ref_pseudos = {
        "C": (psp_dir / f"C.BFD.{file_ext}").resolve(),
        "H": (psp_dir / f"H.BFD.{file_ext}").resolve(),
        "O": (psp_dir / f"O.BFD.{file_ext}").resolve(),
        }

    return (psp_dir, pseudo_list, ref_pseudos)


def test_pp_elem_label():
    from ..pseudopotential import pp_elem_label

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


@isolate_nexus_core
def test_pseudopotentials():
    from ..pseudopotential import Pseudopotentials
    from ..pseudopotential import PseudoFile
    from ..pseudopotential import gamessPPFile

    # empty initialization
    Pseudopotentials()
    PseudoFile()
    gamessPPFile()

    # standard initialization
    file_paths = list(TEST_FILES.values())
    pps = Pseudopotentials(file_paths)

    for fn in TEST_FILES.keys():
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
    from ..developer import obj
    from ..pseudopotential import ppset

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



def test_pseudopotential_classes(tmp_path):
    import numpy as np
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

    qo = qpp.to_obj()
    co = cpp.to_obj()
    del qo.rmin
    del qo.rmax
    assert(object_eq(co,qo,atol=1e-12))

    qmcpack_from_casino_file = tmp_path / 'C.BFD.qmcpack_from_casino'
    cpp.write_qmcpack(qmcpack_from_casino_file)

    qpp_casino = QmcpackPP(qmcpack_from_casino_file)
    assert(object_eq(qpp_casino,qpp,int_as_float=True,atol=1e-12))

#end def test_pseudopotential_classes


def test_read_upf_z_valence(tmp_path):
    upf_v1_file = TEST_FILES["C.BFD.upf"]

    upf_v1_file_zvalence = read_upf_z_valence(upf_v1_file)
    assert(isinstance(upf_v1_file_zvalence, int))
    assert(upf_v1_file_zvalence == 4)

    upf_v1_alt_text = """
<PP_HEADER>
   0         Version Number
   Cu        Element
   NC        Norm - Conserving pseudopotential
    F      Nonlinear Core Correction
SLA  PZ   NOGX NOGC    PZ   Exchange-Correlation functional
 19          Z valence
 0          Total energy
 0.000000   0.000000     Suggested cutoff for wfc and rho
 2           Max angular momentum component
 1163           Number of points in mesh
 3  2     Number of Wavefuncitons, Number of Projectors
 Wavefunctions         nl  l   occ
                       S  0  2.000000
                       P  1  6.000000
                       D  2  10.000000
</PP_HEADER>
"""
    upf_v1_alt_file = tmp_path / "pseudo_alt_v1.upf"
    upf_v1_alt_file.write_text(upf_v1_alt_text)

    upf_v1_alt_file_zvalence = read_upf_z_valence(upf_v1_alt_file)

    assert(isinstance(upf_v1_alt_file_zvalence, int))
    assert(upf_v1_alt_file_zvalence == 19)

    upf_v_201_text = """
<UPF version="2.0.1">
  <PP_INFO>
    ...
    <PP_INPUTFILE>
    ...
    </PP_INPUTFILE>
  </PP_INFO>
  <!--                               -->
  <!-- END OF HUMAN READABLE SECTION -->
  <!--                               -->
    <PP_HEADER
       generated="By Someone"
       author="Someone"
       date="999999"
       comment=""
       element="C "
       pseudo_type="NC"
       relativistic="scalar"
       is_ultrasoft="F"
       is_paw="F"
       is_coulomb="F"
       has_so="F"
       has_wfc="F"
       has_gipaw="F"
       core_correction="F"
       functional="PBE"
       z_valence="    4.00"
       total_psenergy="  0.000E+00"
       rho_cutoff="   6.00000000000E+00"
       l_max="1"
       l_local="-1"
       mesh_size="   8"
       number_of_wfc="0"
       number_of_proj="1"/>
 <PP_MESH>
   <PP_R type="real"  size=" 8" columns="8">
    0.0000    0.0100    0.0200    0.0300    0.0400    0.0500    0.0600    0.0700
   </PP_R>
 </PP_MESH>
"""
    upf_v201_file = tmp_path / "pseudo_v2.upf"
    upf_v201_file.write_text(upf_v_201_text)

    upf_v201_file_zvalence = read_upf_z_valence(upf_v201_file)

    assert(isinstance(upf_v201_file_zvalence, int))
    assert(upf_v201_file_zvalence == 4)

    # Write Z valence with a float return
    upf_v_201_text_float = """
<UPF version="2.0.1">
  <PP_INFO>
    ...
    <PP_INPUTFILE>
    ...
    </PP_INPUTFILE>
  </PP_INFO>
  <!--                               -->
  <!-- END OF HUMAN READABLE SECTION -->
  <!--                               -->
    <PP_HEADER
       generated="By Someone"
       author="Someone"
       date="999999"
       comment=""
       element="C "
       pseudo_type="NC"
       relativistic="scalar"
       is_ultrasoft="F"
       is_paw="F"
       is_coulomb="F"
       has_so="F"
       has_wfc="F"
       has_gipaw="F"
       core_correction="F"
       functional="PBE"
       z_valence="    4.50"
       total_psenergy="  0.000E+00"
       rho_cutoff="   6.00000000000E+00"
       l_max="1"
       l_local="-1"
       mesh_size="   8"
       number_of_wfc="0"
       number_of_proj="1"/>
 <PP_MESH>
   <PP_R type="real"  size=" 8" columns="8">
    0.0000    0.0100    0.0200    0.0300    0.0400    0.0500    0.0600    0.0700
   </PP_R>
 </PP_MESH>
"""
    upf_v201_file_float = tmp_path / "pseudo_float.upf"
    upf_v201_file_float.write_text(upf_v_201_text_float)

    upf_v201_file_float_zvalence = read_upf_z_valence(upf_v201_file_float)

    assert(upf_v201_file_float_zvalence == 4.5)

    # Single-line header
    upf_v_201_text_one_line = """
<PP_HEADER generated="Generated using &quot;atomic&quot; code by A. Dal Corso  v.6.2.2" author="ADC" date=" 2May2018" comment="" element="Zn" pseudo_type="PAW" relativistic="scalar" is_ultrasoft="true" is_paw="true" is_coulomb="false" has_so="false" has_wfc="true" has_gipaw="true" paw_as_gipaw="true" core_correction="true" functional=" SLA  PW   PBX  PBC" z_valence="1.200000000000e1" total_psenergy="-2.434243516297e2" wfc_cutoff="4.363174091908e1" rho_cutoff="2.755329390766e2" l_max="2" l_max_rho="4" l_local="-1" mesh_size="1201" number_of_wfc="3" number_of_proj="6"/>
"""
    upf_v201_file_one_line = tmp_path / "pseudo_one_line.upf"
    upf_v201_file_one_line.write_text(upf_v_201_text_one_line)

    upf_v201_file_one_line_zvalence = read_upf_z_valence(upf_v201_file_one_line)

    assert(upf_v201_file_one_line_zvalence == 12)

    # Strangely formatted header
    upf_v_201_text_strange = """
<PP_HEADER generated="Generated using &quot;atomic&quot; code by A. Dal Corso  v.6.2.2" author="ADC" 

date=" 2May2018" comment="" element="Zn" pseudo_type="PAW" relativistic="scalar" is_ultrasoft="true" is_paw="true" is_coulomb="false" has_so="false" has_wfc="true" has_gipaw="true" 

paw_as_gipaw="true" core_correction="true" functional=" SLA  PW   PBX  PBC" z_valence = "1.200000000000e1" total_psenergy="-2.434243516297e2" 

wfc_cutoff="4.363174091908e1" rho_cutoff="2.755329390766e2" l_max="2" l_max_rho="4" l_local="-1" mesh_size="1201" number_of_wfc="3" number_of_proj="6"/>
"""
    upf_v201_file_strange = tmp_path / "pseudo_strange.upf"
    upf_v201_file_strange.write_text(upf_v_201_text_strange)

    upf_v201_file_strange_zvalence = read_upf_z_valence(upf_v201_file_strange)

    assert(upf_v201_file_strange_zvalence == 12)
#end def test_read_upf_z_valence


def test_read_xml_z_valence(tmp_path):
    xml_file = TEST_FILES["C.BFD.xml"]

    z_valence = read_xml_z_valence(xml_file)

    assert(isinstance(z_valence, int))
    assert(z_valence == 4)

    # Modified version of C.BFD.xml with a float for the Z-valence
    xml_with_float = """
<?xml version="1.0" encoding="UTF-8"?>
<pseudo version="0.5">
  <header symbol="C" atomic-number="6" zval="4.5" relativistic="no" 
   polarized="no" creator="ppconvert" flavor="Troullier-Martins" 
   core-corrections="no" xc-functional-type="GGA" 
   xc-functional-parametrization="Perdew-Burke-Ernzerhof"/>
  <grid type="linear" units="bohr" ri="0" rf="1.00000000000000e+01" 
   npts="10001"/>
  <semilocal units="hartree" format="r*V" npots-down="2" npots-up="0" 
  </semilocal>
  <pseudowave-functions units="electrons/bohr^(-3/2)" 
  </pseudowave-functions>
</pseudo>
"""

    xml_file_float = tmp_path / "pseudo.xml"
    xml_file_float.write_text(xml_with_float)

    z_valence_float = read_xml_z_valence(xml_file_float)

    assert(z_valence_float == 4.5)    
#end def test_read_xml_z_valence


def test_read_potcar_z_valence(tmp_path):
    """POTCAR examples assembled from publicly available information [1]_.

    First test is for a properly formatted POTCAR with an integer-value
    Z-valence. Second test is for an improperly formatted POTCAR with an
    integer-value Z-valence, to test the fallback method. Third test is
    a properly formatted POTCAR with a non-integer-value Z-valence.

    References
    ----------
    .. [1] https://vasp.at/wiki/POTCAR#File_format
    """
    potcar_proper = """
PAW_PBE C 07Sep2000
4.000
TITEL = PAW_PBE C 07Sep2000
LEXCH = PE
POMASS = 47.880; ZVAL = 4.000 mass and valenz
ENMAX = 222.335; ENMIN = 166.751 eV
EAUG = 482.848
Atomic configuration
    8 entries
     n  l   j            E        occ.
     1  0  0.50     -4865.3608   2.0000
     2  0  0.50      -533.1368   2.0000
     2  1  1.50      -440.5031   6.0000
     3  0  0.50       -59.3186   2.0000
     3  1  1.50       -35.7012   6.0000
     3  2  2.50        -1.9157   3.0000
     4  0  0.50        -3.7291   1.0000
     4  3  2.50        -1.3606   0.0000
"""

    potcar_proper_file = tmp_path / "POTCAR"
    potcar_proper_file.write_text(potcar_proper)

    z_valence = read_potcar_z_valence(potcar_proper_file)

    assert(isinstance(z_valence, int))
    assert(z_valence == 4)

    # Modified version of POTCAR with incorrect 2nd line, to test the fallback system
    potcar_improper = """
PAW_PBE C 07Sep2000
########NOT THE RIGHT LINE########
TITEL = PAW_PBE C 07Sep2000
LEXCH = PE
POMASS = 47.880; ZVAL = 4.000 mass and valenz
ENMAX = 222.335; ENMIN = 166.751 eV
EAUG = 482.848
Atomic configuration
    8 entries
     n  l   j            E        occ.
     1  0  0.50     -4865.3608   2.0000
     2  0  0.50      -533.1368   2.0000
     2  1  1.50      -440.5031   6.0000
     3  0  0.50       -59.3186   2.0000
     3  1  1.50       -35.7012   6.0000
     3  2  2.50        -1.9157   3.0000
     4  0  0.50        -3.7291   1.0000
     4  3  2.50        -1.3606   0.0000
"""

    potcar_improper_file = tmp_path / "POTCAR2"
    potcar_improper_file.write_text(potcar_improper)

    z_valence = read_potcar_z_valence(potcar_improper_file)

    assert(isinstance(z_valence, int))
    assert(z_valence == 4)

    # Modified version of POTCAR with a non-integer Z-valence
    potcar_float = """
PAW_PBE C 07Sep2000
4.500
TITEL = PAW_PBE C 07Sep2000
LEXCH = PE
POMASS = 47.880; ZVAL = 4.500 mass and valenz
ENMAX = 222.335; ENMIN = 166.751 eV
EAUG = 482.848
Atomic configuration
    8 entries
     n  l   j            E        occ.
     1  0  0.50     -4865.3608   2.0000
     2  0  0.50      -533.1368   2.0000
     2  1  1.50      -440.5031   6.0000
     3  0  0.50       -59.3186   2.0000
     3  1  1.50       -35.7012   6.0000
     3  2  2.50        -1.9157   3.0000
     4  0  0.50        -3.7291   1.0000
     4  3  2.50        -1.3606   0.0000
"""

    potcar_float_file = tmp_path / "POTCAR3"
    potcar_float_file.write_text(potcar_float)

    z_valence = read_potcar_z_valence(potcar_float_file)

    assert(isinstance(z_valence, float))
    assert(z_valence == 4.5)
#end def test_read_potcar_z_valence


def test_pseudoset_check_code_str():
    for code in ("espresso", "gamess", "vasp", "qmcpack"):
        assert(code == PseudoSet._check_code_str(code))
        assert(code == PseudoSet._check_code_str(code.upper())) # All Uppercase
        assert(code == PseudoSet._check_code_str(code.title())) # Uppercase first letter

    with pytest.warns(
        NexusUserWarning,
        match="Automatically switching code 'pwscf' to 'espresso'"
        ):
        PseudoSet._check_code_str("pwscf")

    with pytest.raises(
        ValueError,
        match="Code 'bean' is not known by Nexus!"
        ):
        PseudoSet._check_code_str("bean")
#end def test_pseudoset_check_code_str


def test_pseudoset_dict(tmp_path):
    xml_dir, _, ref_xml_pseudos = setup_psps(test_dir=tmp_path, code="qmcpack")
    upf_dir, _, ref_upf_pseudos = setup_psps(test_dir=tmp_path, code="espresso")
    gms_dir, _, ref_gms_pseudos = setup_psps(test_dir=tmp_path, code="gamess")
    potcar_dir, _, ref_potcar_pseudos = setup_psps(test_dir=tmp_path, code="vasp")

    xml_pseudo_dict = {
        "C1": ref_xml_pseudos["C"],
        "H1": ref_xml_pseudos["H"],
        "O1": ref_xml_pseudos["O"],
        }
    xml_pseudoset = PseudoSet(
        pseudos = xml_pseudo_dict,
        code    = "qmcpack",
        )

    assert(xml_pseudoset.pseudos     == xml_pseudo_dict)
    assert(xml_pseudoset.pseudo_dirs == {xml_dir})
    assert(xml_pseudoset.code        == "qmcpack")

    upf_pseudo_dict = {
        "C1": ref_upf_pseudos["C"],
        "H1": ref_upf_pseudos["H"],
        "O1": ref_upf_pseudos["O"],
        }
    upf_pseudoset = PseudoSet(
        pseudos = upf_pseudo_dict,
        code    = "espresso",
        )

    assert(upf_pseudoset.pseudos     == upf_pseudo_dict)
    assert(upf_pseudoset.pseudo_dirs == {upf_dir})
    assert(upf_pseudoset.code        == "espresso")

    gms_pseudo_dict = {
        "C1": ref_gms_pseudos["C"],
        "H1": ref_gms_pseudos["H"],
        "O1": ref_gms_pseudos["O"],
        }
    gms_pseudoset = PseudoSet(
        pseudos = gms_pseudo_dict,
        code    = "gamess",
        )

    assert(gms_pseudoset.pseudos     == gms_pseudo_dict)
    assert(gms_pseudoset.pseudo_dirs == {gms_dir})
    assert(gms_pseudoset.code        == "gamess")

    potcar_pseudo_dict = {
        "C1": ref_potcar_pseudos["C"],
        "H1": ref_potcar_pseudos["H"],
        "O1": ref_potcar_pseudos["O"],
        }
    potcar_pseudoset = PseudoSet(
        pseudos = potcar_pseudo_dict,
        code    = "vasp",
        )

    assert(potcar_pseudoset.pseudos     == potcar_pseudo_dict)
    assert(potcar_pseudoset.pseudo_dirs == {potcar_dir})
    assert(potcar_pseudoset.code        == "vasp")
#end def test_pseudoset_dict


def test_pseudoset_list(tmp_path):
    xml_dir, xml_pseudo_list, ref_xml_pseudos = setup_psps(test_dir=tmp_path, code="qmcpack")
    upf_dir, upf_pseudo_list, ref_upf_pseudos = setup_psps(test_dir=tmp_path, code="espresso")
    gms_dir, gms_pseudo_list, ref_gms_pseudos = setup_psps(test_dir=tmp_path, code="gamess")
    potcar_dir, potcar_pseudo_list, ref_potcar_pseudos = setup_psps(test_dir=tmp_path, code="vasp")

    xml_pseudoset = PseudoSet(
        pseudos = xml_pseudo_list,
        code    = "qmcpack",
        )

    assert(xml_pseudoset.pseudos     == ref_xml_pseudos)
    assert(xml_pseudoset.pseudo_dirs == {xml_dir})
    assert(xml_pseudoset.code        == "qmcpack")

    upf_pseudoset = PseudoSet(
        pseudos = upf_pseudo_list,
        code    = "espresso",
        )

    assert(upf_pseudoset.pseudos     == ref_upf_pseudos)
    assert(upf_pseudoset.pseudo_dirs == {upf_dir})
    assert(upf_pseudoset.code        == "espresso")

    gms_pseudoset = PseudoSet(
        pseudos = gms_pseudo_list,
        code    = "gamess",
        )

    assert(gms_pseudoset.pseudos     == ref_gms_pseudos)
    assert(gms_pseudoset.pseudo_dirs == {gms_dir})
    assert(gms_pseudoset.code        == "gamess")

    potcar_pseudoset = PseudoSet(
        pseudos = potcar_pseudo_list,
        code    = "vasp",
        )

    assert(potcar_pseudoset.pseudos     == ref_potcar_pseudos)
    assert(potcar_pseudoset.pseudo_dirs == {potcar_dir})
    assert(potcar_pseudoset.code        == "vasp")
#end def test_pseudoset_list


def test_pseudoset_detect(tmp_path):

    detect_qmcpack = PseudoSet._detect_pseudo_code(["C.BFD.xml", "H.BFD.xml", "O.BFD.XML"])
    assert(detect_qmcpack == "qmcpack")

    detect_espresso = PseudoSet._detect_pseudo_code(["C.BFD.upf", "H.BFD.RRKJ3", "O.BFD.ncpp"])
    assert(detect_espresso == "espresso")

    detect_gamess = PseudoSet._detect_pseudo_code(["C.BFD.gms", "H.BFD.gms", "O.BFD.gms"])
    assert(detect_gamess == "gamess")

    detect_vasp = PseudoSet._detect_pseudo_code(["C.BFD.potcar", "H.BFD.potcar", "O.BFD.potcar"])
    assert(detect_vasp == "vasp")

    with pytest.raises(
        RuntimeError,
        match="Can not detect code from pseudopotential extensions!"
        ):
        PseudoSet._detect_pseudo_code(["C.BFD.xml", "H.BFD.upf", "O.BFD.gms"])

    with pytest.raises(
        RuntimeError,
        match="Can not detect code from pseudopotential extensions!"
        ):
        PseudoSet._detect_pseudo_code(["C.BFD.txt"])

    _, xml_pseudo_list, _ = setup_psps(test_dir=tmp_path, code="qmcpack")
    xml_pseudoset = PseudoSet(pseudos=xml_pseudo_list, code="detect")

    assert(xml_pseudoset.code == "qmcpack")
#end def test_pseudoset_detect


def test_pseudoset_from_dir(tmp_path):
    xml_dir, _, ref_xml_pseudos = setup_psps(test_dir=tmp_path, code="qmcpack")
    upf_dir, _, ref_upf_pseudos = setup_psps(test_dir=tmp_path, code="espresso")
    gms_dir, _, ref_gms_pseudos = setup_psps(test_dir=tmp_path, code="gamess")
    potcar_dir, _, ref_potcar_pseudos = setup_psps(test_dir=tmp_path, code="vasp")

    xml_pseudoset = PseudoSet.from_dir(pseudo_dir=xml_dir, code="qmcpack")

    assert(xml_pseudoset.pseudos     == ref_xml_pseudos)
    assert(xml_pseudoset.pseudo_dirs == {xml_dir})
    assert(xml_pseudoset.code        == "qmcpack")

    upf_pseudoset = PseudoSet.from_dir(pseudo_dir=upf_dir, code="espresso")

    assert(upf_pseudoset.pseudos     == ref_upf_pseudos)
    assert(upf_pseudoset.pseudo_dirs == {upf_dir})
    assert(upf_pseudoset.code        == "espresso")

    gms_pseudoset = PseudoSet.from_dir(pseudo_dir=gms_dir, code="gamess")

    assert(gms_pseudoset.pseudos     == ref_gms_pseudos)
    assert(gms_pseudoset.pseudo_dirs == {gms_dir})
    assert(gms_pseudoset.code        == "gamess")

    potcar_pseudoset = PseudoSet.from_dir(pseudo_dir=potcar_dir, code="vasp")

    assert(potcar_pseudoset.pseudos     == ref_potcar_pseudos)
    assert(potcar_pseudoset.pseudo_dirs == {potcar_dir})
    assert(potcar_pseudoset.code        == "vasp")
#end def test_pseudoset_from_dir


def test_pseudoset_from_dir_detect(tmp_path):
    xml_dir, _, ref_xml_pseudos = setup_psps(test_dir=tmp_path, code="qmcpack")
    upf_dir, _, ref_upf_pseudos = setup_psps(test_dir=tmp_path, code="espresso")
    gms_dir, _, ref_gms_pseudos = setup_psps(test_dir=tmp_path, code="gamess")
    potcar_dir, _, ref_potcar_pseudos = setup_psps(test_dir=tmp_path, code="vasp")

    xml_pseudoset = PseudoSet.from_dir(pseudo_dir=xml_dir, code="detect")

    assert(xml_pseudoset.pseudos     == ref_xml_pseudos)
    assert(xml_pseudoset.pseudo_dirs == {xml_dir})
    assert(xml_pseudoset.code        == "qmcpack")

    upf_pseudoset = PseudoSet.from_dir(pseudo_dir=upf_dir, code="detect")

    assert(upf_pseudoset.pseudos     == ref_upf_pseudos)
    assert(upf_pseudoset.pseudo_dirs == {upf_dir})
    assert(upf_pseudoset.code        == "espresso")

    gms_pseudoset = PseudoSet.from_dir(pseudo_dir=gms_dir, code="detect")

    assert(gms_pseudoset.pseudos     == ref_gms_pseudos)
    assert(gms_pseudoset.pseudo_dirs == {gms_dir})
    assert(gms_pseudoset.code        == "gamess")

    potcar_pseudoset = PseudoSet.from_dir(pseudo_dir=potcar_dir, code="detect")

    assert(potcar_pseudoset.pseudos     == ref_potcar_pseudos)
    assert(potcar_pseudoset.pseudo_dirs == {potcar_dir})
    assert(potcar_pseudoset.code        == "vasp")
#end def test_pseudoset_from_dir_detect


def test_pseudoset_from_dir_filter(tmp_path):
    pseudo_names = (
        "C.BFD.xml",
        "C.BFD.upf",
        "H.BFD.xml",
        "H.BFD.upf",
        "O.BFD.xml",
        "O.BFD.upf",
        "Fe.BFD.ncpp",
        )

    psp_dir = tmp_path / "mixed_pseudos"
    psp_dir.mkdir()
    assert psp_dir.exists(), "Failed to create pseudo directory!"

    extra_dir = psp_dir / "other_directory"
    extra_dir.mkdir()
    assert extra_dir.exists(), "Failed to create extra directory!"

    pseudo_list = []
    for psp in pseudo_names:
        pseudo = (psp_dir / psp).resolve()
        pseudo.touch()
        assert pseudo.exists(), "Failed to create pseudo file!"
        pseudo_list.append(pseudo)

    ref_xml_pseudos = {
        "C": (psp_dir / "C.BFD.xml").resolve(),
        "H": (psp_dir / "H.BFD.xml").resolve(),
        "O": (psp_dir / "O.BFD.xml").resolve(),
        }

    ref_upf_pseudos = {
        "C" : (psp_dir / "C.BFD.upf").resolve(),
        "H" : (psp_dir / "H.BFD.upf").resolve(),
        "O" : (psp_dir / "O.BFD.upf").resolve(),
        "Fe": (psp_dir / "Fe.BFD.ncpp").resolve(),
        }

    xml_pseudoset = PseudoSet.from_dir(
        pseudo_dir  = psp_dir,
        code        = "qmcpack",
        filter_exts = True,
        )

    assert(xml_pseudoset.pseudos     == ref_xml_pseudos)
    assert(xml_pseudoset.pseudo_dirs == {psp_dir})
    assert(xml_pseudoset.code        == "qmcpack")

    upf_pseudoset = PseudoSet.from_dir(
        pseudo_dir  = psp_dir,
        code        = "espresso",
        filter_exts = True,
        )

    assert(upf_pseudoset.pseudos     == ref_upf_pseudos)
    assert(upf_pseudoset.pseudo_dirs == {psp_dir})
    assert(upf_pseudoset.code        == "espresso")

    custom_filter_pseudoset = PseudoSet.from_dir(
        pseudo_dir  = psp_dir,
        code        = "detect",
        filter_exts = [".upf", ".ncpp"],
        )

    assert(custom_filter_pseudoset.pseudos     == ref_upf_pseudos)
    assert(custom_filter_pseudoset.pseudo_dirs == {psp_dir})
    assert(custom_filter_pseudoset.code        == "espresso")
#end def test_pseudoset_from_dir_filter


def test_get_pseudos(tmp_path):
    xml_dir, _, ref_xml_pseudos = setup_psps(test_dir=tmp_path, code="qmcpack")
    xml_pseudoset = PseudoSet.from_dir(pseudo_dir=xml_dir, code="detect")

    elements = ["C", "C", "H"]

    pseudos = xml_pseudoset.get_pseudos(system=elements, code="qmcpack")

    ref_pseudos = {
        "C": ref_xml_pseudos["C"],
        "H": ref_xml_pseudos["H"],
    }
    assert(pseudos == ref_pseudos)

    with pytest.raises(ValueError, match="Tried to get pseudopotentials for"):
        xml_pseudoset.get_pseudos(system=elements, code="espresso")

    with pytest.raises(ValueError, match="No pseudopotential found for label"):
        xml_pseudoset.get_pseudos(system=["C", "H", "Fe"], code="qmcpack")
#end def test_get_pseudos


def test_get_Zeffs():
    psp_dir = TEST_FILES["C.BFD.xml"].parent.resolve()
    ref_Zeffs_default = {
        "C": 4,
        "H": 1,
        "O": 6,
        }

    xml_pseudoset = PseudoSet.from_dir(
        pseudo_dir  = psp_dir,
        code        = "qmcpack",
        filter_exts = True,
        )
    Zeffs = xml_pseudoset.get_Zeffs(elem_labels=["C", "H", "O"])

    assert(Zeffs == ref_Zeffs_default)

    ref_Zeffs_customized = {
        "C": 6,
        "H": 1,
        "O": 6,
        }

    xml_pseudoset_custom = PseudoSet.from_dir(
        pseudo_dir  = psp_dir,
        code        = "qmcpack",
        Zeffs       = {"C": 6},
        filter_exts = True,
        )
    Zeffs = xml_pseudoset_custom.get_Zeffs(elem_labels=["C", "H", "O"])

    assert(Zeffs == ref_Zeffs_customized)

    ref_Zeffs_mixed_ae = {
        "C":   4,
        "H":   1,
        "O":   6,
        "Fe": 26,
        }

    xml_pseudoset_mixed_ae = PseudoSet.from_dir(
        pseudo_dir  = psp_dir,
        code        = "qmcpack",
        filter_exts = True,
        )
    Zeffs = xml_pseudoset_mixed_ae.get_Zeffs(
        elem_labels   = ["C", "H", "O", "Fe"],
        missing_as_ae = True,
        )

    assert(Zeffs == ref_Zeffs_mixed_ae)

    with pytest.raises(ValueError, match="No pseudopotential found for label"):
        xml_pseudoset.get_Zeffs(elem_labels=["C", "H", "Fe"])

    with pytest.raises(ValueError, match="Can not determine element for label"):
        xml_pseudoset.get_Zeffs(elem_labels=["C", "H", "NotAnElement"], missing_as_ae=True)
#end def test_get_Zeffs
