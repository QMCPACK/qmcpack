import pytest
from copy import deepcopy
from . import NexusTestOrder
pytestmark = pytest.mark.order(NexusTestOrder.PSEUDOPOTENTIAL)

from ..generic import generic_settings, NexusError, NexusUserWarning
generic_settings.raise_error = True

import numpy as np
from pathlib import Path
from typing import Literal
from . import isolate_nexus_core, TEST_DIR
from ..testing import value_eq,object_eq
from nexus.pseudo_set import read_upf_z_valence, read_qmcpack_xml_z_valence, read_potcar_z_valence
from nexus.pseudo_set import PseudoSet, ppset
from nexus.nexus_base import nexus_core
from nexus.physical_system import generate_physical_system


TEST_FILES = {
    "C.BFD.gms": TEST_DIR / "test_pseudopotential_files/C.BFD.gms",
    "C.BFD.upf": TEST_DIR / "../examples/qmcpack/pseudopotentials/C.BFD.upf",
    "C.BFD.xml": TEST_DIR / "../examples/qmcpack/pseudopotentials/C.BFD.xml",
    }

for file in TEST_FILES.values():
    assert(file.exists()), f"Test file not found! {file}"


def setup_psps(
    test_dir: Path,
    code: Literal["espresso", "gamess", "vasp", "qmcpack", "rmg", "pyscf"],
    ) -> tuple[Path, list[Path], dict[str, Path]]:
    """Take a test's temp directory and populate with dummy pseudopotential files."""
    match code:
        # At most 3 extensions can be added
        case "espresso":
            file_ext = (".upf", ".ncpp")
        case "qmcpack":
            file_ext = ".xml"
        case "gamess":
            file_ext = ".gms"
        case "vasp":
            file_ext = "POTCAR"
        case "rmg":
            file_ext = (".upf", ".xml")
        case "pyscf":
            file_ext = (".nwchem", ".gth")
        case _:
            msg = (
                "Invalid call to `setup_for_pseudoset()`!\n"
                f"Code supplied is {code}, but must be one of: {', '.join(PseudoSet.known_codes)}"
                )
            raise pytest.UsageError(msg)
    if file_ext == "POTCAR":
        pseudo_names = (
            "C/POTCAR",
            "H/POTCAR",
            "O/POTCAR",
            )
    elif isinstance(file_ext, tuple):
        pseudo_names = (
            f"C.BFD.{file_ext[0]}",
            f"H.BFD.{file_ext[1]}",
            f"O.BFD.{file_ext[-1]}",
            )
    else:
        pseudo_names = (
            f"C.BFD.{file_ext}",
            f"H.BFD.{file_ext}",
            f"O.BFD.{file_ext}",
            )

    psp_dir = test_dir / f"{code}_pseudos"
    psp_dir.mkdir()
    assert psp_dir.exists(), "Failed to create pseudo directory!"

    pseudo_list = []
    for psp in pseudo_names:
        if "POTCAR" in psp:
            potcar_dir = psp_dir / psp.split("/")[0]
            potcar_dir.mkdir()
            assert potcar_dir.exists(), "Failed to create POTCAR directory!"

        pseudo = (psp_dir / psp).resolve()
        pseudo.touch()
        assert pseudo.exists(), "Failed to create pseudo file!"
        pseudo_list.append(pseudo)

    if file_ext == "POTCAR":
        ref_pseudos = {
            "C": (psp_dir / "C" / "POTCAR").resolve(),
            "H": (psp_dir / "H" / "POTCAR").resolve(),
            "O": (psp_dir / "O" / "POTCAR").resolve(),
            }
    elif isinstance(file_ext, tuple):
        ref_pseudos = {
            "C": (psp_dir / f"C.BFD.{file_ext[0]}").resolve(),
            "H": (psp_dir / f"H.BFD.{file_ext[1]}").resolve(),
            "O": (psp_dir / f"O.BFD.{file_ext[-1]}").resolve(),
            }
    else:
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


def test_read_qmcpack_xml_z_valence(tmp_path):
    xml_file = TEST_FILES["C.BFD.xml"]

    z_valence = read_qmcpack_xml_z_valence(xml_file)

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

    z_valence_float = read_qmcpack_xml_z_valence(xml_file_float)

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
    for code in ("espresso", "gamess", "vasp", "qmcpack", "rmg", "pyscf"):
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
    qmcpack_dir,  _, ref_qmcpack_pseudos  = setup_psps(test_dir=tmp_path, code="qmcpack")
    espresso_dir, _, ref_espresso_pseudos = setup_psps(test_dir=tmp_path, code="espresso")
    gamess_dir,   _, ref_gamess_pseudos   = setup_psps(test_dir=tmp_path, code="gamess")
    vasp_dir,     _, ref_vasp_pseudos     = setup_psps(test_dir=tmp_path, code="vasp")
    rmg_dir,      _, ref_rmg_pseudos      = setup_psps(test_dir=tmp_path, code="rmg")
    pyscf_dir,    _, ref_pyscf_pseudos    = setup_psps(test_dir=tmp_path, code="pyscf")

    qmcpack_pseudo_dict = {
        "C1": ref_qmcpack_pseudos["C"],
        "H1": ref_qmcpack_pseudos["H"],
        "O1": ref_qmcpack_pseudos["O"],
        }
    qmcpack_pseudoset = PseudoSet(
        pseudos = qmcpack_pseudo_dict,
        codes   = "qmcpack",
        )

    assert(qmcpack_pseudoset.pseudos     == qmcpack_pseudo_dict)
    assert(qmcpack_pseudoset.pseudo_dirs == {qmcpack_dir})
    assert(qmcpack_pseudoset.codes       == {"qmcpack"})

    espresso_pseudo_dict = {
        "C1": ref_espresso_pseudos["C"],
        "H1": ref_espresso_pseudos["H"],
        "O1": ref_espresso_pseudos["O"],
        }
    espresso_pseudoset = PseudoSet(
        pseudos = espresso_pseudo_dict,
        codes   = "espresso",
        )

    assert(espresso_pseudoset.pseudos     == espresso_pseudo_dict)
    assert(espresso_pseudoset.pseudo_dirs == {espresso_dir})
    assert(espresso_pseudoset.codes       == {"espresso"})

    gamess_pseudo_dict = {
        "C1": ref_gamess_pseudos["C"],
        "H1": ref_gamess_pseudos["H"],
        "O1": ref_gamess_pseudos["O"],
        }
    gamess_pseudoset = PseudoSet(
        pseudos = gamess_pseudo_dict,
        codes   = "gamess",
        )

    assert(gamess_pseudoset.pseudos     == gamess_pseudo_dict)
    assert(gamess_pseudoset.pseudo_dirs == {gamess_dir})
    assert(gamess_pseudoset.codes       == {"gamess"})

    vasp_pseudo_dict = {
        "C": ref_vasp_pseudos["C"],
        "H": ref_vasp_pseudos["H"],
        "O": ref_vasp_pseudos["O"],
        }
    vasp_pseudoset = PseudoSet(
        pseudos = vasp_pseudo_dict,
        codes   = "vasp",
        )

    assert(vasp_pseudoset.pseudos     == vasp_pseudo_dict)
    assert(vasp_pseudoset.pseudo_dirs == {vasp_dir})
    assert(vasp_pseudoset.codes       == {"vasp"})

    rmg_pseudo_dict = {
        "C": ref_rmg_pseudos["C"],
        "H": ref_rmg_pseudos["H"],
        "O": ref_rmg_pseudos["O"],
        }
    rmg_pseudoset = PseudoSet(
        pseudos = rmg_pseudo_dict,
        codes   = "rmg",
        )

    assert(rmg_pseudoset.pseudos     == rmg_pseudo_dict)
    assert(rmg_pseudoset.pseudo_dirs == {rmg_dir})
    assert(rmg_pseudoset.codes       == {"rmg"})

    pyscf_pseudo_dict = {
        "C": ref_pyscf_pseudos["C"],
        "H": ref_pyscf_pseudos["H"],
        "O": ref_pyscf_pseudos["O"],
        }
    pyscf_pseudoset = PseudoSet(
        pseudos = pyscf_pseudo_dict,
        codes   = "pyscf",
        )

    assert(pyscf_pseudoset.pseudos     == pyscf_pseudo_dict)
    assert(pyscf_pseudoset.pseudo_dirs == {pyscf_dir})
    assert(pyscf_pseudoset.codes       == {"pyscf"})
#end def test_pseudoset_dict


def test_pseudoset_list(tmp_path):
    qmcpack_dir,  qmcpack_pseudo_list,  ref_qmcpack_pseudos  = setup_psps(test_dir=tmp_path, code="qmcpack")
    espresso_dir, espresso_pseudo_list, ref_espresso_pseudos = setup_psps(test_dir=tmp_path, code="espresso")
    gamess_dir,   gamess_pseudo_list,   ref_gamess_pseudos   = setup_psps(test_dir=tmp_path, code="gamess")
    vasp_dir,     vasp_pseudo_list,     ref_vasp_pseudos     = setup_psps(test_dir=tmp_path, code="vasp")
    rmg_dir,      rmg_pseudo_list,      ref_rmg_pseudos      = setup_psps(test_dir=tmp_path, code="rmg")
    pyscf_dir,    pyscf_pseudo_list,    ref_pyscf_pseudos    = setup_psps(test_dir=tmp_path, code="pyscf")

    qmcpack_pseudoset = PseudoSet(
        pseudos = qmcpack_pseudo_list,
        codes    = "qmcpack",
        )

    assert(qmcpack_pseudoset.pseudos     == ref_qmcpack_pseudos)
    assert(qmcpack_pseudoset.pseudo_dirs == {qmcpack_dir})
    assert(qmcpack_pseudoset.codes       == {"qmcpack"})

    espresso_pseudoset = PseudoSet(
        pseudos = espresso_pseudo_list,
        codes    = "espresso",
        )

    assert(espresso_pseudoset.pseudos     == ref_espresso_pseudos)
    assert(espresso_pseudoset.pseudo_dirs == {espresso_dir})
    assert(espresso_pseudoset.codes       == {"espresso"})

    gamess_pseudoset = PseudoSet(
        pseudos = gamess_pseudo_list,
        codes    = "gamess",
        )

    assert(gamess_pseudoset.pseudos     == ref_gamess_pseudos)
    assert(gamess_pseudoset.pseudo_dirs == {gamess_dir})
    assert(gamess_pseudoset.codes       == {"gamess"})

    vasp_pseudoset = PseudoSet(
        pseudos = vasp_pseudo_list,
        codes    = "vasp",
        )

    assert(vasp_pseudoset.pseudos     == ref_vasp_pseudos)
    assert(vasp_pseudoset.pseudo_dirs == {vasp_dir})
    assert(vasp_pseudoset.codes       == {"vasp"})

    rmg_pseudoset = PseudoSet(
        pseudos = rmg_pseudo_list,
        codes    = "rmg",
        )

    assert(rmg_pseudoset.pseudos     == ref_rmg_pseudos)
    assert(rmg_pseudoset.pseudo_dirs == {rmg_dir})
    assert(rmg_pseudoset.codes       == {"rmg"})

    pyscf_pseudoset = PseudoSet(
        pseudos = pyscf_pseudo_list,
        codes   = "pyscf",
        )

    assert(pyscf_pseudoset.pseudos     == ref_pyscf_pseudos)
    assert(pyscf_pseudoset.pseudo_dirs == {pyscf_dir})
    assert(pyscf_pseudoset.codes       == {"pyscf"})
#end def test_pseudoset_list


def test_pseudoset_detect(tmp_path):

    detect_qmcpack_rmg = PseudoSet._detect_pseudo_code(["C.BFD.xml", "H.BFD.xml", "O.BFD.XML"])
    assert(detect_qmcpack_rmg == {"qmcpack", "rmg"})

    detect_espresso = PseudoSet._detect_pseudo_code(["C.BFD.upf", "H.BFD.RRKJ3", "O.BFD.ncpp"])
    assert(detect_espresso == {"espresso"})

    detect_espresso_rmg = PseudoSet._detect_pseudo_code(["C.BFD.upf", "H.BFD.upf", "O.BFD.upf"])
    assert(detect_espresso_rmg == {"espresso", "rmg"})

    detect_gamess = PseudoSet._detect_pseudo_code(["C.BFD.gms", "H.BFD.gms", "O.BFD.gamess"])
    assert(detect_gamess == {"gamess"})

    detect_vasp = PseudoSet._detect_pseudo_code(["C/POTCAR", "H/POTCAR", "O.vasp"])
    assert(detect_vasp == {"vasp"})

    detect_pyscf = PseudoSet._detect_pseudo_code(["C.BFD.nwchem", "H.BFD.nwchem", "O.BFD.gth"])
    assert(detect_pyscf == {"pyscf"})

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
    qmcpack_pseudoset = PseudoSet(pseudos=xml_pseudo_list, codes="detect")

    assert(qmcpack_pseudoset.codes == {"qmcpack", "rmg"})
#end def test_pseudoset_detect


def test_pseudoset_from_dir(tmp_path):
    qmcpack_dir,  _, ref_qmcpack_pseudos  = setup_psps(test_dir=tmp_path, code="qmcpack")
    espresso_dir, _, ref_espresso_pseudos = setup_psps(test_dir=tmp_path, code="espresso")
    gamess_dir,   _, ref_gamess_pseudos   = setup_psps(test_dir=tmp_path, code="gamess")
    vasp_dir,     _, ref_vasp_pseudos     = setup_psps(test_dir=tmp_path, code="vasp")
    rmg_dir,      _, ref_rmg_pseudos      = setup_psps(test_dir=tmp_path, code="rmg")
    pyscf_dir,    _, ref_pyscf_pseudos    = setup_psps(test_dir=tmp_path, code="pyscf")

    qmcpack_pseudoset = PseudoSet.from_dir(pseudo_dir=qmcpack_dir, code ="qmcpack")

    assert(qmcpack_pseudoset.pseudos     == ref_qmcpack_pseudos)
    assert(qmcpack_pseudoset.pseudo_dirs == {qmcpack_dir})
    assert(qmcpack_pseudoset.codes       == {"qmcpack"})

    espresso_pseudoset = PseudoSet.from_dir(pseudo_dir=espresso_dir, code ="espresso")

    assert(espresso_pseudoset.pseudos     == ref_espresso_pseudos)
    assert(espresso_pseudoset.pseudo_dirs == {espresso_dir})
    assert(espresso_pseudoset.codes       == {"espresso"})

    gamess_pseudoset = PseudoSet.from_dir(pseudo_dir=gamess_dir, code ="gamess")

    assert(gamess_pseudoset.pseudos     == ref_gamess_pseudos)
    assert(gamess_pseudoset.pseudo_dirs == {gamess_dir})
    assert(gamess_pseudoset.codes       == {"gamess"})

    vasp_pseudoset = PseudoSet.from_dir(pseudo_dir=vasp_dir, code ="vasp")

    assert(vasp_pseudoset.pseudos     == ref_vasp_pseudos)
    assert(vasp_pseudoset.pseudo_dirs == {vasp_dir})
    assert(vasp_pseudoset.codes       == {"vasp"})

    rmg_pseudoset = PseudoSet.from_dir(pseudo_dir=rmg_dir, code ="rmg")

    assert(rmg_pseudoset.pseudos     == ref_rmg_pseudos)
    assert(rmg_pseudoset.pseudo_dirs == {rmg_dir})
    assert(rmg_pseudoset.codes       == {"rmg"})

    pyscf_pseudoset = PseudoSet.from_dir(pseudo_dir=pyscf_dir, code ="pyscf")

    assert(pyscf_pseudoset.pseudos     == ref_pyscf_pseudos)
    assert(pyscf_pseudoset.pseudo_dirs == {pyscf_dir})
    assert(pyscf_pseudoset.codes       == {"pyscf"})
#end def test_pseudoset_from_dir


def test_pseudoset_from_dir_detect(tmp_path):
    qmcpack_dir,  _, ref_qmcpack_pseudos  = setup_psps(test_dir=tmp_path, code="qmcpack")
    espresso_dir, _, ref_espresso_pseudos = setup_psps(test_dir=tmp_path, code="espresso")
    gamess_dir,   _, ref_gamess_pseudos   = setup_psps(test_dir=tmp_path, code="gamess")
    vasp_dir,     _, ref_vasp_pseudos     = setup_psps(test_dir=tmp_path, code="vasp")
    rmg_dir,      _, ref_rmg_pseudos      = setup_psps(test_dir=tmp_path, code="rmg")
    pyscf_dir,    _, ref_pyscf_pseudos    = setup_psps(test_dir=tmp_path, code="pyscf")

    qmcpack_pseudoset = PseudoSet.from_dir(pseudo_dir=qmcpack_dir, code ="detect")

    assert(qmcpack_pseudoset.pseudos     == ref_qmcpack_pseudos)
    assert(qmcpack_pseudoset.pseudo_dirs == {qmcpack_dir})
    assert(qmcpack_pseudoset.codes       == {"qmcpack", "rmg"})

    espresso_pseudoset = PseudoSet.from_dir(pseudo_dir=espresso_dir, code ="detect")

    assert(espresso_pseudoset.pseudos     == ref_espresso_pseudos)
    assert(espresso_pseudoset.pseudo_dirs == {espresso_dir})
    assert(espresso_pseudoset.codes       == {"espresso"})

    gamess_pseudoset = PseudoSet.from_dir(pseudo_dir=gamess_dir, code ="detect")

    assert(gamess_pseudoset.pseudos     == ref_gamess_pseudos)
    assert(gamess_pseudoset.pseudo_dirs == {gamess_dir})
    assert(gamess_pseudoset.codes       == {"gamess"})

    vasp_pseudoset = PseudoSet.from_dir(pseudo_dir=vasp_dir, code ="detect")

    assert(vasp_pseudoset.pseudos     == ref_vasp_pseudos)
    assert(vasp_pseudoset.pseudo_dirs == {vasp_dir})
    assert(vasp_pseudoset.codes       == {"vasp"})

    rmg_pseudoset = PseudoSet.from_dir(pseudo_dir=rmg_dir, code ="detect")

    assert(rmg_pseudoset.pseudos     == ref_rmg_pseudos)
    assert(rmg_pseudoset.pseudo_dirs == {rmg_dir})
    assert(rmg_pseudoset.codes       == {"rmg"})

    pyscf_pseudoset = PseudoSet.from_dir(pseudo_dir=pyscf_dir, code ="detect")

    assert(pyscf_pseudoset.pseudos     == ref_pyscf_pseudos)
    assert(pyscf_pseudoset.pseudo_dirs == {pyscf_dir})
    assert(pyscf_pseudoset.codes       == {"pyscf"})
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

    ref_qmcpack_pseudos = {
        "C": (psp_dir / "C.BFD.xml").resolve(),
        "H": (psp_dir / "H.BFD.xml").resolve(),
        "O": (psp_dir / "O.BFD.xml").resolve(),
        }

    ref_espresso_pseudos = {
        "C" : (psp_dir / "C.BFD.upf").resolve(),
        "H" : (psp_dir / "H.BFD.upf").resolve(),
        "O" : (psp_dir / "O.BFD.upf").resolve(),
        "Fe": (psp_dir / "Fe.BFD.ncpp").resolve(),
        }

    qmcpack_pseudoset = PseudoSet.from_dir(
        pseudo_dir  = psp_dir,
        code        = "qmcpack",
        ext_filter  = None,
        )

    assert(qmcpack_pseudoset.pseudos     == ref_qmcpack_pseudos)
    assert(qmcpack_pseudoset.pseudo_dirs == {psp_dir})
    assert(qmcpack_pseudoset.codes       == {"qmcpack"})

    espresso_pseudoset = PseudoSet.from_dir(
        pseudo_dir  = psp_dir,
        code        = "espresso",
        ext_filter  = None,
        )

    assert(espresso_pseudoset.pseudos     == ref_espresso_pseudos)
    assert(espresso_pseudoset.pseudo_dirs == {psp_dir})
    assert(espresso_pseudoset.codes       == {"espresso"})

    custom_filter_pseudoset = PseudoSet.from_dir(
        pseudo_dir  = psp_dir,
        code        = "detect",
        ext_filter  = [".upf", ".ncpp"],
        )

    assert(custom_filter_pseudoset.pseudos     == ref_espresso_pseudos)
    assert(custom_filter_pseudoset.pseudo_dirs == {psp_dir})
    assert(custom_filter_pseudoset.codes       == {"espresso"})
#end def test_pseudoset_from_dir_filter


def test_pseudoset_from_dir_pattern(tmp_path):
    pseudo_names = (
        "C.BFD.xml",
        "C.BFD.upf",
        "H.BFD.xml",
        "H.BFD.upf",
        "O.BFD.xml",
        "O.BFD.upf",
        "Fe.BFD.ncpp",
        "C_ONCV_PBE-1.2.upf",
        "H_ONCV_PBE-1.2.upf",
        "O_ONCV_PBE-1.2.upf",
        "Fe_ONCV_PBE-1.2.upf",
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

    ref_qmcpack_pseudos = {
        "C": (psp_dir / "C.BFD.xml").resolve(),
        "H": (psp_dir / "H.BFD.xml").resolve(),
        "O": (psp_dir / "O.BFD.xml").resolve(),
        }

    ref_espresso_bfd_pseudos = {
        "C" : (psp_dir / "C.BFD.upf").resolve(),
        "H" : (psp_dir / "H.BFD.upf").resolve(),
        "O" : (psp_dir / "O.BFD.upf").resolve(),
        "Fe": (psp_dir / "Fe.BFD.ncpp").resolve(),
        }

    ref_espresso_oncv_pseudos = {
        "C" : (psp_dir / "C_ONCV_PBE-1.2.upf").resolve(),
        "H" : (psp_dir / "H_ONCV_PBE-1.2.upf").resolve(),
        "O" : (psp_dir / "O_ONCV_PBE-1.2.upf").resolve(),
        "Fe": (psp_dir / "Fe_ONCV_PBE-1.2.upf").resolve(),
        }

    qmcpack_pseudoset = PseudoSet.from_dir(
        pseudo_dir  = psp_dir,
        code        = "qmcpack",
        ext_filter  = None,
        )

    assert(qmcpack_pseudoset.pseudos     == ref_qmcpack_pseudos)
    assert(qmcpack_pseudoset.pseudo_dirs == {psp_dir})
    assert(qmcpack_pseudoset.codes       == {"qmcpack"})

    espresso_bfd_pseudoset = PseudoSet.from_dir(
        pseudo_dir  = psp_dir,
        code        = "espresso",
        ext_filter  = None,
        pattern     = r"\.BFD\.",
        )

    assert(espresso_bfd_pseudoset.pseudos     == ref_espresso_bfd_pseudos)
    assert(espresso_bfd_pseudoset.pseudo_dirs == {psp_dir})
    assert(espresso_bfd_pseudoset.codes       == {"espresso"})

    espresso_oncv_pseudoset = PseudoSet.from_dir(
        pseudo_dir  = psp_dir,
        code        = "espresso",
        ext_filter  = None,
        pattern     = r"_ONCV_PBE-1\.2",
        )

    assert(espresso_oncv_pseudoset.pseudos     == ref_espresso_oncv_pseudos)
    assert(espresso_oncv_pseudoset.pseudo_dirs == {psp_dir})
    assert(espresso_oncv_pseudoset.codes       == {"espresso"})
#end def test_pseudoset_from_dir_pattern


def test_vasp_pattern_exclude(tmp_path):
    pseudo_names = (
        "C/POTCAR",
        "H/POTCAR",
        "O/POTCAR",
        "C_sv/POTCAR",
        "H_sv/POTCAR",
        "O_sv/POTCAR",
        "C_GW/POTCAR",
        "H_GW/POTCAR",
        "O_GW/POTCAR",
        "C_sv_GW/POTCAR",
        "H_sv_GW/POTCAR",
        "O_sv_GW/POTCAR",
        )

    psp_dir = tmp_path / "vasp_pseudos"
    psp_dir.mkdir()
    assert psp_dir.exists(), "Failed to create pseudo directory!"

    pseudo_list = []
    for psp in pseudo_names:
        potcar_dir = psp_dir / psp.split("/")[0]
        potcar_dir.mkdir()
        assert potcar_dir.exists(), "Failed to create POTCAR directory!"

        pseudo = (psp_dir / psp).resolve()
        pseudo.touch()
        assert pseudo.exists(), "Failed to create pseudo file!"
        pseudo_list.append(pseudo)

    ref_reg_pseudos = {
        "C": (psp_dir / "C" / "POTCAR").resolve(),
        "H": (psp_dir / "H" / "POTCAR").resolve(),
        "O": (psp_dir / "O" / "POTCAR").resolve(),
    }
    ref_sv_pseudos = {
        "C": (psp_dir / "C_sv" / "POTCAR").resolve(),
        "H": (psp_dir / "H_sv" / "POTCAR").resolve(),
        "O": (psp_dir / "O_sv" / "POTCAR").resolve(),
    }
    ref_gw_pseudos = {
        "C": (psp_dir / "C_GW" / "POTCAR").resolve(),
        "H": (psp_dir / "H_GW" / "POTCAR").resolve(),
        "O": (psp_dir / "O_GW" / "POTCAR").resolve(),
    }
    ref_sv_gw_pseudos = {
        "C": (psp_dir / "C_sv_GW" / "POTCAR").resolve(),
        "H": (psp_dir / "H_sv_GW" / "POTCAR").resolve(),
        "O": (psp_dir / "O_sv_GW" / "POTCAR").resolve(),
    }

    reg_pseudoset = PseudoSet.from_dir(
        pseudo_dir = psp_dir,
        code       = "vasp",
        pattern    = r"^((?!_).)*$",
    )

    assert(reg_pseudoset.pseudos == ref_reg_pseudos)

    sv_pseudoset = PseudoSet.from_dir(
        pseudo_dir = psp_dir,
        code       = "vasp",
        pattern    = r"_sv$",
    )

    assert(sv_pseudoset.pseudos == ref_sv_pseudos)

    gw_pseudoset = PseudoSet.from_dir(
        pseudo_dir = psp_dir,
        code       = "vasp",
        pattern    = r"(?<!sv)_GW",
    )

    assert(gw_pseudoset.pseudos == ref_gw_pseudos)

    sv_gw_pseudoset = PseudoSet.from_dir(
        pseudo_dir = psp_dir,
        code       = "vasp",
        pattern    = r"_sv_GW",
    )

    assert(sv_gw_pseudoset.pseudos == ref_sv_gw_pseudos)
#end def test_vasp_pattern_exclude


def test_pseudoset_from_mixed_dir(tmp_path):
    pseudo_names = (
        "C.BFD.xml",
        "H.BFD.xml",
        "O.BFD.xml",
        "C.BFD.ncpp",
        "H.BFD.ncpp",
        "O.BFD.ncpp",
        "Fe.BFD.upf",
        "C.BFD.gms",
        "H.BFD.gms",
        "O.BFD.gamess",
        "C.BFD.nwchem",
        "H.BFD.nwchem",
        "O.BFD.gth",
        "C/POTCAR",
        "H/POTCAR",
        "O/POTCAR",
        )

    psp_dir = tmp_path / "mixed_pseudos"
    psp_dir.mkdir()
    assert psp_dir.exists(), "Failed to create pseudo directory!"

    pseudo_list = []
    for psp in pseudo_names:
        if "POTCAR" in psp:
            potcar_dir = psp_dir / psp.split("/")[0]
            potcar_dir.mkdir()
            assert potcar_dir.exists(), "Failed to create POTCAR directory!"

        pseudo = (psp_dir / psp).resolve()
        pseudo.touch()
        assert pseudo.exists(), "Failed to create pseudo file!"
        pseudo_list.append(pseudo)

    ref_pseudos = {
        "qmcpack": PseudoSet(
            pseudos = {
                "C": (psp_dir / "C.BFD.xml").resolve(),
                "H": (psp_dir / "H.BFD.xml").resolve(),
                "O": (psp_dir / "O.BFD.xml").resolve(),
                },
            codes = {"qmcpack"},
            ),
        "espresso": PseudoSet(
            pseudos = {
                "C" : (psp_dir / "C.BFD.ncpp").resolve(),
                "H" : (psp_dir / "H.BFD.ncpp").resolve(),
                "O" : (psp_dir / "O.BFD.ncpp").resolve(),
                "Fe": (psp_dir / "Fe.BFD.upf").resolve(),
                },
            codes={"espresso"},
            ),
        "gamess": PseudoSet(
            pseudos = {
                "C": (psp_dir / "C.BFD.gms").resolve(),
                "H": (psp_dir / "H.BFD.gms").resolve(),
                "O": (psp_dir / "O.BFD.gamess").resolve(),
                },
            codes={"gamess"},
            ),
        "vasp": PseudoSet(
            pseudos = {
                "C": (psp_dir / "C" / "POTCAR").resolve(),
                "H": (psp_dir / "H" / "POTCAR").resolve(),
                "O": (psp_dir / "O" / "POTCAR").resolve(),
                },
            codes={"vasp"},
            ),
        "rmg": PseudoSet(
            pseudos = {
                "C": (psp_dir / "C.BFD.xml").resolve(),
                "H": (psp_dir / "H.BFD.xml").resolve(),
                "O": (psp_dir / "O.BFD.xml").resolve(),
                "Fe": (psp_dir / "Fe.BFD.upf").resolve(),
                },
            codes={"rmg"},
            ),
        "pyscf": PseudoSet(
            pseudos = {
                "C": (psp_dir / "C.BFD.nwchem").resolve(),
                "H": (psp_dir / "H.BFD.nwchem").resolve(),
                "O": (psp_dir / "O.BFD.gth").resolve(),
                },
            codes={"pyscf"},
            ),
        }

    pseudoset = PseudoSet.from_mixed_dir(
        pseudo_dir = psp_dir,
        codes      = None, # Default values
        extensions = None, # Default values
        )

    assert(pseudoset.keys() == ref_pseudos.keys())

    qmcpack_calc = pseudoset["qmcpack"]
    qmcpack_ref = ref_pseudos["qmcpack"]
    assert(qmcpack_calc.pseudos     == qmcpack_ref.pseudos)
    assert(qmcpack_calc.codes       == qmcpack_ref.codes)
    assert(qmcpack_calc.pseudo_dirs == qmcpack_ref.pseudo_dirs)

    espresso_calc = pseudoset["espresso"]
    espresso_ref = ref_pseudos["espresso"]
    assert(espresso_calc.pseudos     == espresso_ref.pseudos)
    assert(espresso_calc.codes       == espresso_ref.codes)
    assert(espresso_calc.pseudo_dirs == espresso_ref.pseudo_dirs)

    gamess_calc = pseudoset["gamess"]
    gamess_ref = ref_pseudos["gamess"]
    assert(gamess_calc.pseudos     == gamess_ref.pseudos)
    assert(gamess_calc.codes       == gamess_ref.codes)
    assert(gamess_calc.pseudo_dirs == gamess_ref.pseudo_dirs)

    vasp_calc = pseudoset["vasp"]
    vasp_ref = ref_pseudos["vasp"]
    assert(vasp_calc.pseudos     == vasp_ref.pseudos)
    assert(vasp_calc.codes       == vasp_ref.codes)
    assert(vasp_calc.pseudo_dirs == vasp_ref.pseudo_dirs)

    rmg_calc = pseudoset["rmg"]
    rmg_ref = ref_pseudos["rmg"]
    assert(rmg_calc.pseudos     == rmg_ref.pseudos)
    assert(rmg_calc.codes       == rmg_ref.codes)
    assert(rmg_calc.pseudo_dirs == rmg_ref.pseudo_dirs)

    pyscf_calc = pseudoset["pyscf"]
    pyscf_ref = ref_pseudos["pyscf"]
    assert(pyscf_calc.pseudos     == pyscf_ref.pseudos)
    assert(pyscf_calc.codes       == pyscf_ref.codes)
    assert(pyscf_calc.pseudo_dirs == pyscf_ref.pseudo_dirs)

    pseudoset = PseudoSet.from_mixed_dir(
        pseudo_dir = psp_dir,
        codes      = ["qmcpack", "espresso"],
        extensions = None,
        )

    assert(set(pseudoset.keys()) == {"qmcpack", "espresso"})

    qmcpack_calc = pseudoset["qmcpack"]
    qmcpack_ref = ref_pseudos["qmcpack"]
    assert(qmcpack_calc.pseudos     == qmcpack_ref.pseudos)
    assert(qmcpack_calc.codes        == qmcpack_ref.codes)
    assert(qmcpack_calc.pseudo_dirs == qmcpack_ref.pseudo_dirs)

    espresso_calc = pseudoset["espresso"]
    espresso_ref = ref_pseudos["espresso"]
    assert(espresso_calc.pseudos     == espresso_ref.pseudos)
    assert(espresso_calc.codes        == espresso_ref.codes)
    assert(espresso_calc.pseudo_dirs == espresso_ref.pseudo_dirs)

    pseudoset = PseudoSet.from_mixed_dir(
        pseudo_dir = psp_dir,
        codes      = {"qmcpack", "espresso"},
        extensions = {
            "qmcpack": {".xml"},
            "espresso": {".ncpp"},
            },
        )

    assert(set(pseudoset.keys()) == {"espresso", "qmcpack"})

    qmcpack_calc = pseudoset["qmcpack"]
    qmcpack_ref = ref_pseudos["qmcpack"]
    assert(qmcpack_calc.pseudos     == qmcpack_ref.pseudos)
    assert(qmcpack_calc.codes       == qmcpack_ref.codes)
    assert(qmcpack_calc.pseudo_dirs == qmcpack_ref.pseudo_dirs)

    espresso_calc = pseudoset["espresso"]
    espresso_ref = ref_pseudos["espresso"]

    del espresso_ref.pseudos["Fe"] # Only grabbed `.ncpp` files

    assert(espresso_calc.pseudos     == espresso_ref.pseudos)
    assert(espresso_calc.codes       == espresso_ref.codes)
    assert(espresso_calc.pseudo_dirs == espresso_ref.pseudo_dirs)
#end def test_pseudoset_from_mixed_dir


def test_pseudoset_from_mixed_dir_rmg_collision(tmp_path):
    pseudo_names = (
        "C.BFD.xml",
        "H.BFD.xml",
        "O.BFD.xml",
        "C.BFD.upf",
        "H.BFD.upf",
        "O.BFD.upf",
        # Throw in some misc files to make sure they get ignored
        "C.tmp",
        "C2.other",
        )

    psp_dir = tmp_path / "mixed_pseudos"
    psp_dir.mkdir()
    assert psp_dir.exists(), "Failed to create pseudo directory!"

    pseudo_list = []
    for psp in pseudo_names:
        pseudo = (psp_dir / psp).resolve()
        pseudo.touch()
        assert pseudo.exists(), "Failed to create pseudo file!"
        pseudo_list.append(pseudo)

    ref_pseudos = {
        "qmcpack": PseudoSet(
            pseudos = {
                "C": (psp_dir / "C.BFD.xml").resolve(),
                "H": (psp_dir / "H.BFD.xml").resolve(),
                "O": (psp_dir / "O.BFD.xml").resolve(),
                },
            codes = {"qmcpack"},
            ),
        "espresso": PseudoSet(
            pseudos = {
                "C" : (psp_dir / "C.BFD.upf").resolve(),
                "H" : (psp_dir / "H.BFD.upf").resolve(),
                "O" : (psp_dir / "O.BFD.upf").resolve(),
                },
            codes={"espresso"},
            ),
        "gamess": PseudoSet(
            pseudos = {},
            codes={"gamess"},
            ),
        "vasp": PseudoSet(
            pseudos = {},
            codes={"vasp"},
            ),
        "rmg": PseudoSet(
            pseudos = {
                "C" : (psp_dir / "C.BFD.upf").resolve(),
                "H" : (psp_dir / "H.BFD.upf").resolve(),
                "O" : (psp_dir / "O.BFD.upf").resolve(),
                },
            codes={"rmg"},
            ),
        "pyscf": PseudoSet(
            pseudos = {},
            codes={"pyscf"},
            ),
        }

    with pytest.raises(
        RuntimeError,
        match=(
            "Duplicate element detected for code 'rmg'\n"
            "Either remove 'rmg' from the selected codes, or specify "
            "`filters` and/or `patterns` to ensure the collision does not happen"
            ),
        ):
        pseudoset = PseudoSet.from_mixed_dir(
            pseudo_dir = psp_dir,
            codes      = None, # Default values
            extensions = None, # Default values
            )

    # Make sure we can follow the instructions in the message from the error.
    pseudoset = PseudoSet.from_mixed_dir(
        pseudo_dir = psp_dir,
        codes      = ["qmcpack", "espresso"],
        extensions = None,
        )

    assert(set(pseudoset.keys()) == {"qmcpack", "espresso"})

    qmcpack_calc = pseudoset["qmcpack"]
    qmcpack_ref = ref_pseudos["qmcpack"]
    assert(qmcpack_calc.pseudos     == qmcpack_ref.pseudos)
    assert(qmcpack_calc.codes       == qmcpack_ref.codes)
    assert(qmcpack_calc.pseudo_dirs == qmcpack_ref.pseudo_dirs)

    espresso_calc = pseudoset["espresso"]
    espresso_ref = ref_pseudos["espresso"]
    assert(espresso_calc.pseudos     == espresso_ref.pseudos)
    assert(espresso_calc.codes       == espresso_ref.codes)
    assert(espresso_calc.pseudo_dirs == espresso_ref.pseudo_dirs)

    # Filter for rmg, leave codes as None
    pseudoset = PseudoSet.from_mixed_dir(
        pseudo_dir = psp_dir,
        codes      = None,
        extensions = {
            "rmg": {".upf"}
            },
        )

    assert(set(pseudoset.keys()) == {"espresso", "gamess", "vasp", "qmcpack", "rmg", "pyscf"})

    qmcpack_calc = pseudoset["qmcpack"]
    qmcpack_ref = ref_pseudos["qmcpack"]
    assert(qmcpack_calc.pseudos     == qmcpack_ref.pseudos)
    assert(qmcpack_calc.codes       == qmcpack_ref.codes)
    assert(qmcpack_calc.pseudo_dirs == qmcpack_ref.pseudo_dirs)

    espresso_calc = pseudoset["espresso"]
    espresso_ref = ref_pseudos["espresso"]
    assert(espresso_calc.pseudos     == espresso_ref.pseudos)
    assert(espresso_calc.codes       == espresso_ref.codes)
    assert(espresso_calc.pseudo_dirs == espresso_ref.pseudo_dirs)

    gamess_calc = pseudoset["gamess"]
    gamess_ref = ref_pseudos["gamess"]
    assert(gamess_calc.pseudos     == gamess_ref.pseudos)
    assert(gamess_calc.codes       == gamess_ref.codes)
    assert(gamess_calc.pseudo_dirs == gamess_ref.pseudo_dirs)

    vasp_calc = pseudoset["vasp"]
    vasp_ref = ref_pseudos["vasp"]
    assert(vasp_calc.pseudos     == vasp_ref.pseudos)
    assert(vasp_calc.codes       == vasp_ref.codes)
    assert(vasp_calc.pseudo_dirs == vasp_ref.pseudo_dirs)

    rmg_calc = pseudoset["rmg"]
    rmg_ref = ref_pseudos["rmg"]
    assert(rmg_calc.pseudos     == rmg_ref.pseudos)
    assert(rmg_calc.codes       == rmg_ref.codes)
    assert(rmg_calc.pseudo_dirs == rmg_ref.pseudo_dirs)

    pyscf_calc = pseudoset["pyscf"]
    pyscf_ref = ref_pseudos["pyscf"]
    assert(pyscf_calc.pseudos     == pyscf_ref.pseudos)
    assert(pyscf_calc.codes       == pyscf_ref.codes)
    assert(pyscf_calc.pseudo_dirs == pyscf_ref.pseudo_dirs)

    # Pattern for rmg, leave codes as None
    pseudoset = PseudoSet.from_mixed_dir(
        pseudo_dir = psp_dir,
        codes      = None,
        patterns   = {
            "rmg": "upf",
            },
        )

    assert(set(pseudoset.keys()) == {"espresso", "gamess", "vasp", "qmcpack", "rmg", "pyscf"})

    qmcpack_calc = pseudoset["qmcpack"]
    qmcpack_ref = ref_pseudos["qmcpack"]
    assert(qmcpack_calc.pseudos     == qmcpack_ref.pseudos)
    assert(qmcpack_calc.codes       == qmcpack_ref.codes)
    assert(qmcpack_calc.pseudo_dirs == qmcpack_ref.pseudo_dirs)

    espresso_calc = pseudoset["espresso"]
    espresso_ref = ref_pseudos["espresso"]
    assert(espresso_calc.pseudos     == espresso_ref.pseudos)
    assert(espresso_calc.codes       == espresso_ref.codes)
    assert(espresso_calc.pseudo_dirs == espresso_ref.pseudo_dirs)

    gamess_calc = pseudoset["gamess"]
    gamess_ref = ref_pseudos["gamess"]
    assert(gamess_calc.pseudos     == gamess_ref.pseudos)
    assert(gamess_calc.codes       == gamess_ref.codes)
    assert(gamess_calc.pseudo_dirs == gamess_ref.pseudo_dirs)

    vasp_calc = pseudoset["vasp"]
    vasp_ref = ref_pseudos["vasp"]
    assert(vasp_calc.pseudos     == vasp_ref.pseudos)
    assert(vasp_calc.codes       == vasp_ref.codes)
    assert(vasp_calc.pseudo_dirs == vasp_ref.pseudo_dirs)

    rmg_calc = pseudoset["rmg"]
    rmg_ref = ref_pseudos["rmg"]
    assert(rmg_calc.pseudos     == rmg_ref.pseudos)
    assert(rmg_calc.codes       == rmg_ref.codes)
    assert(rmg_calc.pseudo_dirs == rmg_ref.pseudo_dirs)

    pyscf_calc = pseudoset["pyscf"]
    pyscf_ref = ref_pseudos["pyscf"]
    assert(pyscf_calc.pseudos     == pyscf_ref.pseudos)
    assert(pyscf_calc.codes       == pyscf_ref.codes)
    assert(pyscf_calc.pseudo_dirs == pyscf_ref.pseudo_dirs)
#end def test_pseudoset_from_mixed_dir_rmg_collision


def test_pseudoset_from_mixed_dir_espresso_collision(tmp_path):
    pseudo_names = (
        "C.BFD.upf",
        "C1.uspp.UPF",
        "C_4.ncpp",
        )

    psp_dir = tmp_path / "espresso_collision"
    psp_dir.mkdir()
    assert psp_dir.exists(), "Failed to create pseudo directory!"

    pseudo_list = []
    for psp in pseudo_names:
        pseudo = (psp_dir / psp).resolve()
        pseudo.touch()
        assert pseudo.exists(), "Failed to create pseudo file!"
        pseudo_list.append(pseudo)

    with pytest.raises(
        RuntimeError,
        match=(
            "Duplicate element detected for code 'espresso'\n"
            "Either remove 'espresso' from the selected codes, or specify "
            "`filters` and/or `patterns` to ensure the collision does not happen"
            ),
        ):
        _ = PseudoSet.from_mixed_dir(
            pseudo_dir = psp_dir,
            codes      = None, # Default values
            extensions = None, # Default values
            )

    # Filter by pattern
    _ = PseudoSet.from_mixed_dir(
        pseudo_dir = psp_dir,
        codes      = {"espresso"},
        patterns   = {"espresso": "BFD"},
        )

    _ = PseudoSet.from_mixed_dir(
        pseudo_dir = psp_dir,
        codes      = {"espresso"},
        patterns   = {"espresso": "uspp"},
        )

    # Filter by extension
    _ = PseudoSet.from_mixed_dir(
        pseudo_dir = psp_dir,
        codes      = {"espresso"},
        extensions = {"espresso": ".ncpp"},
        )
#end def test_pseudoset_from_mixed_dir_espresso_collision


def test_pseudoset_from_mixed_dir_vasp_collision(tmp_path):
    pseudo_names = (
        "C/POTCAR",
        "N/POTCAR",
        "C4-special.vasp",
        )

    psp_dir = tmp_path / "vasp_collision"
    psp_dir.mkdir()
    assert psp_dir.exists(), "Failed to create pseudo directory!"

    pseudo_list = []
    for psp in pseudo_names:
        if "POTCAR" in psp:
            potcar_dir = psp_dir / psp.split("/")[0]
            potcar_dir.mkdir()
            assert potcar_dir.exists(), "Failed to create POTCAR directory!"

        pseudo = (psp_dir / psp).resolve()
        pseudo.touch()
        assert pseudo.exists(), "Failed to create pseudo file!"
        pseudo_list.append(pseudo)

    with pytest.raises(
            RuntimeError,
            match=(
                "Duplicate element detected for code 'vasp'\n"
                "Either remove 'vasp' from the selected codes, or specify "
                "`filters` and/or `patterns` to ensure the collision does not happen"
                ),
            ):
            _ = PseudoSet.from_mixed_dir(
                pseudo_dir = psp_dir,
                codes      = {"vasp"},
                extensions = None,
                )
    
    # Filter by pattern
    _ = PseudoSet.from_mixed_dir(
        pseudo_dir = psp_dir,
        codes      = {"vasp"},
        patterns   = {"vasp": "special"},
        )

    # Filter by extension
    _ = PseudoSet.from_mixed_dir(
        pseudo_dir = psp_dir,
        codes      = {"vasp"},
        extensions = {"vasp": ".vasp"},
        )

    _ = PseudoSet.from_mixed_dir(
        pseudo_dir = psp_dir,
        codes      = {"vasp"},
        extensions = {"vasp": "POTCAR"},
        )

    _ = PseudoSet.from_mixed_dir(
        pseudo_dir = psp_dir,
        codes      = {"vasp"},
        extensions = {"vasp": "potcar"},
        )
#end def test_pseudoset_from_mixed_dir_vasp_collision


def test_pseudoset_from_mixed_dir_gamess_collision(tmp_path):
    pseudo_names = (
        "C.gms",
        "C_special.gamess",
        "N.gamess",
        )

    psp_dir = tmp_path / "gamess_collision"
    psp_dir.mkdir()
    assert psp_dir.exists(), "Failed to create pseudo directory!"

    pseudo_list = []
    for psp in pseudo_names:
        pseudo = (psp_dir / psp).resolve()
        pseudo.touch()
        assert pseudo.exists(), "Failed to create pseudo file!"
        pseudo_list.append(pseudo)

    with pytest.raises(
            RuntimeError,
            match=(
                "Duplicate element detected for code 'gamess'\n"
                "Either remove 'gamess' from the selected codes, or specify "
                "`filters` and/or `patterns` to ensure the collision does not happen"
                ),
            ):
            _ = PseudoSet.from_mixed_dir(
                pseudo_dir = psp_dir,
                codes      = None, # Default values
                extensions = None, # Default values
                )

    with pytest.raises(
            RuntimeError,
            match=(
                "Duplicate element detected for code 'gamess'\n"
                "Either remove 'gamess' from the selected codes, or specify "
                "`filters` and/or `patterns` to ensure the collision does not happen"
                ),
            ):
            _ = PseudoSet.from_mixed_dir(
                pseudo_dir = psp_dir,
                codes      = {"gamess"},
                extensions = None,
                )
    
    # Filter by pattern
    _ = PseudoSet.from_mixed_dir(
        pseudo_dir = psp_dir,
        codes      = {"gamess"},
        patterns   = {"gamess": "special"},
        )

    # Filter by extension
    _ = PseudoSet.from_mixed_dir(
        pseudo_dir = psp_dir,
        codes      = {"gamess"},
        extensions = {"gamess": ".gamess"},
        )

    _ = PseudoSet.from_mixed_dir(
        pseudo_dir = psp_dir,
        codes      = {"gamess"},
        extensions = {"gamess": ".gms"},
        )
#end def test_pseudoset_from_mixed_dir_gamess_collision


def test_pseudoset_from_mixed_dir_pyscf_collision(tmp_path):
    pseudo_names = (
        "C.nwchem",
        "C_special.gth",
        "N.gth",
        )

    psp_dir = tmp_path / "pyscf_collision"
    psp_dir.mkdir()
    assert psp_dir.exists(), "Failed to create pseudo directory!"

    pseudo_list = []
    for psp in pseudo_names:
        pseudo = (psp_dir / psp).resolve()
        pseudo.touch()
        assert pseudo.exists(), "Failed to create pseudo file!"
        pseudo_list.append(pseudo)

    with pytest.raises(
            RuntimeError,
            match=(
                "Duplicate element detected for code 'pyscf'\n"
                "Either remove 'pyscf' from the selected codes, or specify "
                "`filters` and/or `patterns` to ensure the collision does not happen"
                ),
            ):
            _ = PseudoSet.from_mixed_dir(
                pseudo_dir = psp_dir,
                codes      = None, # Default values
                extensions = None, # Default values
                )

    with pytest.raises(
            RuntimeError,
            match=(
                "Duplicate element detected for code 'pyscf'\n"
                "Either remove 'pyscf' from the selected codes, or specify "
                "`filters` and/or `patterns` to ensure the collision does not happen"
                ),
            ):
            _ = PseudoSet.from_mixed_dir(
                pseudo_dir = psp_dir,
                codes      = {"pyscf"},
                extensions = None,
                )
    
    # Filter by pattern
    _ = PseudoSet.from_mixed_dir(
        pseudo_dir = psp_dir,
        codes      = {"pyscf"},
        patterns   = {"pyscf": "special"},
        )

    # Filter by extension
    _ = PseudoSet.from_mixed_dir(
        pseudo_dir = psp_dir,
        codes      = {"pyscf"},
        extensions = {"pyscf": ".gth"},
        )

    _ = PseudoSet.from_mixed_dir(
        pseudo_dir = psp_dir,
        codes      = {"pyscf"},
        extensions = {"pyscf": ".nwchem"},
        )
#end def test_pseudoset_from_mixed_dir_pyscf_collision


def test_pseudoset_from_mixed_dir_mismatches(tmp_path):
    with pytest.raises(
        FileNotFoundError,
        match="Can not find pseudopotential directory:",
        ):
        PseudoSet.from_mixed_dir(pseudo_dir="/path/to/nowhere")

    tmp_file = tmp_path / "random_file.bean"
    tmp_file.touch()

    with pytest.raises(
        NotADirectoryError,
        match="Specified path does not point to a directory:",
        ):
        PseudoSet.from_mixed_dir(pseudo_dir=tmp_file)

    with pytest.raises(
        ValueError,
        match="Mismatch between provided filters and codes!",
        ):
        PseudoSet.from_mixed_dir(
            pseudo_dir = tmp_path,
            codes      = {"espresso"},
            extensions = {
                "espresso": ".upf",
                "rmg": ".upf", # Not in the specified codes
                },
            )

    with pytest.raises(
        ValueError,
        match="Mismatch between provided code Zeff map and codes!",
        ):
        PseudoSet.from_mixed_dir(
            pseudo_dir    = tmp_path,
            codes         = {"espresso"},
            code_Zeff_map = {
                "espresso": {"H": 1},
                "rmg": {"H": 1}, # Not in the specified codes
                },
            )

    with pytest.raises(
        ValueError,
        match="Mismatch between provided patterns and codes!",
        ):
        PseudoSet.from_mixed_dir(
            pseudo_dir = tmp_path,
            codes      = {"espresso"},
            patterns   = {
                "espresso": r"BFD",
                "rmg": r"BFD", # Not in the specified codes
                },
            )
#end def test_pseudoset_from_mixed_dir_mismatches


def test_get_pseudos(tmp_path):
    qmcpack_dir, _, ref_qmcpack_pseudos = setup_psps(test_dir=tmp_path, code="qmcpack")
    qmcpack_pseudoset = PseudoSet.from_dir(pseudo_dir=qmcpack_dir, code ="detect")
    ref_pseudos = {ref_qmcpack_pseudos["C"], ref_qmcpack_pseudos["H"]}

    elements = ["C", "C", "H"]

    pseudos = qmcpack_pseudoset.get_pseudos(system=elements, code="qmcpack")

    assert(pseudos == ref_pseudos)

    system = generate_physical_system(
        elem = ["C", "H", "H", "H", "H"],
        pos = np.empty((5,3), dtype=float),
        C = 4,
        H = 1,
        )

    pseudos = qmcpack_pseudoset.get_pseudos(system=system, code="qmcpack")

    assert(pseudos == ref_pseudos)

    with pytest.raises(ValueError, match="Tried to get pseudopotentials for"):
        qmcpack_pseudoset.get_pseudos(system=elements, code="espresso")

    with pytest.raises(ValueError, match="No pseudopotential found for label"):
        qmcpack_pseudoset.get_pseudos(system=["C", "H", "Fe"], code="qmcpack")
#end def test_get_pseudos


def test_get_Zeff():
    psp_dir = TEST_FILES["C.BFD.xml"].parent.resolve()
    ref_Zeff_default = {
        "C": 4,
        "H": 1,
        "O": 6,
        }

    qmcpack_pseudoset = PseudoSet.from_dir(
        pseudo_dir  = psp_dir,
        code        = "qmcpack",
        ext_filter  = None,
        )
    Zeff = qmcpack_pseudoset.get_Zeff(elem_labels=["C", "H", "O"])

    assert(Zeff == ref_Zeff_default)

    ref_Zeff_customized = {
        "C": 6,
        "H": 1,
        "O": 6,
        }

    qmcpack_pseudoset_custom = PseudoSet.from_dir(
        pseudo_dir  = psp_dir,
        code        = "qmcpack",
        Zeff_map    = {"C": 6},
        ext_filter  = None,
        )
    Zeff = qmcpack_pseudoset_custom.get_Zeff(elem_labels=["C", "H", "O"])

    assert(Zeff == ref_Zeff_customized)

    ref_Zeff_mixed_ae = {
        "C":   4,
        "H":   1,
        "O":   6,
        "Fe": 26,
        }

    qmcpack_pseudoset_mixed_ae = PseudoSet.from_dir(
        pseudo_dir  = psp_dir,
        code        = "qmcpack",
        ext_filter  = None,
        )
    Zeff = qmcpack_pseudoset_mixed_ae.get_Zeff(
        elem_labels   = ["C", "H", "O", "Fe"],
        missing_as_ae = True,
        )

    assert(Zeff == ref_Zeff_mixed_ae)

    with pytest.raises(ValueError, match="No pseudopotential found for label"):
        qmcpack_pseudoset.get_Zeff(elem_labels=["C", "H", "Fe"])

    with pytest.raises(ValueError, match="Can not determine element for label"):
        qmcpack_pseudoset.get_Zeff(elem_labels=["C", "H", "NotAnElement"], missing_as_ae=True)
#end def test_get_Zeff


@isolate_nexus_core
def test_register_legacy_ppset(tmp_path):
    pseudo_names = (
        "C.BFD.xml",
        "H.BFD.xml",
        "O.BFD.xml",
        "C.BFD.upf",
        "H.BFD.upf",
        "O.BFD.upf",
        "C.BFD.gms",
        "H.BFD.gms",
        "O.BFD.gms",
        )

    psp_dir = tmp_path / "mixed_pseudos"
    psp_dir.mkdir()
    assert psp_dir.exists(), "Failed to create pseudo directory!"

    nexus_core.file_locations += [str(psp_dir)]

    pseudo_list = []
    for psp in pseudo_names:
        pseudo = (psp_dir / psp).resolve()
        pseudo.touch()
        assert pseudo.exists(), "Failed to create pseudo file!"
        pseudo_list.append(pseudo)

    pseudo_files = PseudoSet.pseudo_files
    labeled_pseudosets = PseudoSet.labeled_pseudosets
    PseudoSet.pseudo_files = {
        pseudo.name:str(pseudo) for pseudo in pseudo_list
        }
    PseudoSet.labeled_pseudosets = {}

    ref_pseudos = {
        "bfd": {
            "qmcpack": PseudoSet(
                pseudos = {
                    "C": (psp_dir / "C.BFD.xml").resolve(),
                    "H": (psp_dir / "H.BFD.xml").resolve(),
                    "O": (psp_dir / "O.BFD.xml").resolve(),
                    },
                ),
            "espresso": PseudoSet(
                pseudos = {
                    "C" : (psp_dir / "C.BFD.upf").resolve(),
                    "H" : (psp_dir / "H.BFD.upf").resolve(),
                    "O" : (psp_dir / "O.BFD.upf").resolve(),
                    },
                ),
            "gamess": PseudoSet(
                pseudos = {
                    "C": (psp_dir / "C.BFD.gms").resolve(),
                    "H": (psp_dir / "H.BFD.gms").resolve(),
                    "O": (psp_dir / "O.BFD.gms").resolve(),
                    },
                ),
            }
        }

    ppset(
        label   = 'bfd',
        pwscf   = ["C.BFD.upf", "H.BFD.upf", "O.BFD.upf"],
        qmcpack = ["C.BFD.xml", "H.BFD.xml", "O.BFD.xml"],
        gamess  = ["C.BFD.gms", "H.BFD.gms", "O.BFD.gms"],
        )

    assert set(PseudoSet.labeled_pseudosets) == {'bfd'}
    assert set(PseudoSet.labeled_pseudosets['bfd']) == {
        'espresso','gamess','qmcpack'
        }

    for code,ref_pseudoset in ref_pseudos['bfd'].items():
        pseudoset = PseudoSet.labeled_pseudosets['bfd'][code]
        assert pseudoset.pseudos == ref_pseudoset.pseudos
        assert pseudoset.codes == ref_pseudoset.codes
        assert pseudoset.pseudo_dirs == ref_pseudoset.pseudo_dirs
    #end for

    system = generate_physical_system(
        elem = ['O','C','H'],
        pos  = np.empty((3,3),dtype=float),
        C    = 4,
        H    = 1,
        O    = 6,
        )
    remapped = PseudoSet.pseudo_remap('qmcpack','bfd',system)
    assert list(remapped) == ['C.BFD.xml','H.BFD.xml','O.BFD.xml']
    assert remapped == {
        filename:PseudoSet.pseudo_files[filename] for filename in remapped
        }

    explicit = ['O.BFD.xml','C.BFD.xml']
    remapped = PseudoSet.pseudo_remap('qmcpack',explicit,system)
    assert list(remapped) == explicit

    with pytest.raises(NexusError,match='label "bfd" is already registered'):
        ppset(label='bfd',qmcpack=['C.BFD.xml','H.BFD.xml','O.BFD.xml'])

    with pytest.raises(NexusError,match='are not present in PseudoSet.pseudo_files'):
        ppset(label='missing',qmcpack=['Ne.missing.xml'])

    with pytest.raises(NexusError,match='same elements'):
        ppset(
            label   = 'unequal_elements',
            pwscf   = ['C.BFD.upf','H.BFD.upf'],
            qmcpack = ['C.BFD.xml','O.BFD.xml'],
            )

    with pytest.raises(NexusError,match='not compatible with that code'):
        ppset(label='incompatible',gamess=['C.BFD.xml','H.BFD.xml','O.BFD.xml'])

    with pytest.raises(NexusError,match='label "unknown" is not registered'):
        PseudoSet.pseudo_remap('qmcpack','unknown',system)

    with pytest.raises(NexusError,match='are not present in PseudoSet.pseudo_files'):
        PseudoSet.pseudo_remap('qmcpack',['Ne.missing.xml'],system)

    missing_species = generate_physical_system(
        elem = ['C','N'],
        pos  = np.empty((2,3),dtype=float),
        C    = 4,
        N    = 5,
        )
    with pytest.raises(NexusError,match=r"does not contain species \['N'\]"):
        PseudoSet.pseudo_remap('qmcpack','bfd',missing_species)

    ppset(
        label = 'shared_upf',
        pwscf = ['C.BFD.upf','H.BFD.upf','O.BFD.upf'],
        rmg   = ['C.BFD.upf','H.BFD.upf','O.BFD.upf'],
        )
    shared = PseudoSet.labeled_pseudosets['shared_upf']
    assert set(shared) == {'espresso','rmg'}
    assert shared['espresso'] is not shared['rmg']

    PseudoSet.pseudo_files = pseudo_files
    PseudoSet.labeled_pseudosets = labeled_pseudosets
#end def test_register_legacy_ppset


def test_pseudoset_repr(tmp_path):
    qmcpack_dir,  _, ref_qmcpack_pseudos  = setup_psps(test_dir=tmp_path, code="qmcpack")
    pseudoset = PseudoSet(
        pseudos  =  ref_qmcpack_pseudos,
        codes    = "qmcpack",
        Zeff_map = {
            "H": 1,
            "C": 4,
            "O": 6,
            }
        )

    ref_repr = f"""\
PseudoSet(
    codes = {{'qmcpack'}},
    pseudos = {{
        'C': PosixPath('{qmcpack_dir}/C.BFD..xml'),
        'H': PosixPath('{qmcpack_dir}/H.BFD..xml'),
        'O': PosixPath('{qmcpack_dir}/O.BFD..xml'),
    }},
    Zeff_map = {{
        'H': 1,
        'C': 4,
        'O': 6,
    }},
)
"""
    assert(repr(pseudoset) == ref_repr)
#end def test_pseudoset_repr
