import pytest

from . import NexusTestOrder

pytestmark = pytest.mark.order(NexusTestOrder.PSEUDOSET)

import re
from copy import deepcopy
from pathlib import Path
from types import MappingProxyType
from typing import Literal

import numpy as np

from nexus.nexus_base import nexus_core
from nexus.physical_system import generate_physical_system
from nexus.pseudoset import PseudoSet, ppset, generate_pseudoset
from nexus.pseudoset import read_potcar_z_valence, read_qmcpack_xml_z_valence, read_upf_z_valence

from ..generic import NexusUserWarning
from . import TEST_DIR, isolate_nexus_core

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
#end def setup_psps


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
        ('Co' , 'Co' , Path('/path/to/pseudo_dir/Co-ONCV.upf') ),
        ]

    for reflabel,refsymbol,ppfile in ppfiles:
        label,symbol,is_elem = pp_elem_label(ppfile)
        assert(is_elem)
        assert(label==reflabel)
        assert(symbol==refsymbol)
    #end for

    with pytest.raises(RuntimeError,match=r"file name .* is invalid"):
        pp_elem_label('../C.xml',guard=True)

    with pytest.raises(
        RuntimeError,
        match="Cannot determine element for pseudopotential file",
        ):
        pp_elem_label('invalid.xml',guard=True)
#end def test_pp_elem_label


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


def test_from_dir(tmp_path):
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
#end def test_from_dir


def test_from_dir_detect(tmp_path):
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
#end def test_from_dir_detect


def test_from_dir_filter(tmp_path):
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
        pseudo_dir = psp_dir,
        code       = "qmcpack",
        extension  = None,
        )

    assert(qmcpack_pseudoset.pseudos     == ref_qmcpack_pseudos)
    assert(qmcpack_pseudoset.pseudo_dirs == {psp_dir})
    assert(qmcpack_pseudoset.codes       == {"qmcpack"})

    espresso_pseudoset = PseudoSet.from_dir(
        pseudo_dir = psp_dir,
        code       = "espresso",
        extension  = None,
        )

    assert(espresso_pseudoset.pseudos     == ref_espresso_pseudos)
    assert(espresso_pseudoset.pseudo_dirs == {psp_dir})
    assert(espresso_pseudoset.codes       == {"espresso"})

    custom_filter_pseudoset = PseudoSet.from_dir(
        pseudo_dir = psp_dir,
        code       = "detect",
        extension  = [".upf", ".ncpp"],
        )

    assert(custom_filter_pseudoset.pseudos     == ref_espresso_pseudos)
    assert(custom_filter_pseudoset.pseudo_dirs == {psp_dir})
    assert(custom_filter_pseudoset.codes       == {"espresso"})
#end def test_from_dir_filter


def test_from_dir_include(tmp_path):
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
        pseudo_dir = psp_dir,
        code       = "qmcpack",
        extension  = None,
        )

    assert(qmcpack_pseudoset.pseudos     == ref_qmcpack_pseudos)
    assert(qmcpack_pseudoset.pseudo_dirs == {psp_dir})
    assert(qmcpack_pseudoset.codes       == {"qmcpack"})

    espresso_bfd_pseudoset = PseudoSet.from_dir(
        pseudo_dir = psp_dir,
        code       = "espresso",
        extension  = None,
        include    = "*BFD*",
        )

    assert(espresso_bfd_pseudoset.pseudos     == ref_espresso_bfd_pseudos)
    assert(espresso_bfd_pseudoset.pseudo_dirs == {psp_dir})
    assert(espresso_bfd_pseudoset.codes       == {"espresso"})

    espresso_oncv_pseudoset = PseudoSet.from_dir(
        pseudo_dir = psp_dir,
        code       = "espresso",
        extension  = None,
        include    = "*_ONCV_PBE-1.2*",
        )

    assert(espresso_oncv_pseudoset.pseudos     == ref_espresso_oncv_pseudos)
    assert(espresso_oncv_pseudoset.pseudo_dirs == {psp_dir})
    assert(espresso_oncv_pseudoset.codes       == {"espresso"})
#end def test_from_dir_include


def test_vasp_exclude(tmp_path):
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
        exclude    = "*_*",
    )

    assert(reg_pseudoset.pseudos == ref_reg_pseudos)

    sv_pseudoset = PseudoSet.from_dir(
        pseudo_dir = psp_dir,
        code       = "vasp",
        include    = "*_sv",
    )

    assert(sv_pseudoset.pseudos == ref_sv_pseudos)

    gw_pseudoset = PseudoSet.from_dir(
        pseudo_dir = psp_dir,
        code       = "vasp",
        include    = "*_GW",
        exclude    = "*sv*",
    )

    assert(gw_pseudoset.pseudos == ref_gw_pseudos)

    sv_gw_pseudoset = PseudoSet.from_dir(
        pseudo_dir = psp_dir,
        code       = "vasp",
        include    = "*_sv_GW",
    )

    assert(sv_gw_pseudoset.pseudos == ref_sv_gw_pseudos)
#end def test_vasp_exclude


def test_from_mixed_dir(tmp_path):
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
#end def test_from_mixed_dir


def test_from_mixed_dir_rmg_collision(tmp_path):
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
        ValueError,
        match=re.compile(
            pattern=(
                r"Can not provide multiple pseudos for the same element!.*"
                r"Problem detected for code 'rmg'"
                ),
            flags=re.DOTALL,
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
#end def test_from_mixed_dir_rmg_collision


def test_from_mixed_dir_espresso_collision(tmp_path):
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
        ValueError,
        match=re.compile(
            pattern=(
                r"Can not provide multiple pseudos for the same element!.*"
                r"Problem detected for code 'espresso'"
                ),
            flags=re.DOTALL,
            ),
        ):
        _ = PseudoSet.from_mixed_dir(
            pseudo_dir = psp_dir,
            codes      = None, # Default values
            extensions = None, # Default values
            )

    # Use include to only select some
    _ = PseudoSet.from_mixed_dir(
        pseudo_dir = psp_dir,
        codes      = {"espresso"},
        include    = {"espresso": "*BFD*"},
        )

    _ = PseudoSet.from_mixed_dir(
        pseudo_dir = psp_dir,
        codes      = {"espresso"},
        include    = {"espresso": "*uspp*"},
        )

    # Filter by extension
    _ = PseudoSet.from_mixed_dir(
        pseudo_dir = psp_dir,
        codes      = {"espresso"},
        extensions = {"espresso": ".ncpp"},
        )
#end def test_from_mixed_dir_espresso_collision


def test_from_mixed_dir_vasp_collision(tmp_path):
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
        ValueError,
        match=re.compile(
            pattern=(
                r"Can not provide multiple pseudos for the same element!.*"
                r"Problem detected for code 'vasp'"
                ),
            flags=re.DOTALL,
            ),
        ):
            _ = PseudoSet.from_mixed_dir(
                pseudo_dir = psp_dir,
                codes      = {"vasp"},
                extensions = None,
                )

    # Filter by include
    _ = PseudoSet.from_mixed_dir(
        pseudo_dir = psp_dir,
        codes      = {"vasp"},
        include    = {"vasp": "special"},
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
#end def test_from_mixed_dir_vasp_collision


def test_from_mixed_dir_gamess_collision(tmp_path):
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
        ValueError,
        match=re.compile(
            pattern=(
                r"Can not provide multiple pseudos for the same element!.*"
                r"Problem detected for code 'gamess'"
                ),
            flags=re.DOTALL,
            ),
        ):
            _ = PseudoSet.from_mixed_dir(
                pseudo_dir = psp_dir,
                codes      = None, # Default values
                extensions = None, # Default values
                )

    with pytest.raises(
        ValueError,
        match=re.compile(
            pattern=(
                r"Can not provide multiple pseudos for the same element!.*"
                r"Problem detected for code 'gamess'"
                ),
            flags=re.DOTALL,
            ),
        ):
            _ = PseudoSet.from_mixed_dir(
                pseudo_dir = psp_dir,
                codes      = {"gamess"},
                extensions = None,
                )

    # Filter by include
    _ = PseudoSet.from_mixed_dir(
        pseudo_dir = psp_dir,
        codes      = {"gamess"},
        include    = {"gamess": "special"},
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
#end def test_from_mixed_dir_gamess_collision


def test_from_mixed_dir_pyscf_collision(tmp_path):
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
        ValueError,
        match=re.compile(
            pattern=(
                r"Can not provide multiple pseudos for the same element!.*"
                r"Problem detected for code 'pyscf'"
                ),
            flags=re.DOTALL,
            ),
        ):
            _ = PseudoSet.from_mixed_dir(
                pseudo_dir = psp_dir,
                codes      = None, # Default values
                extensions = None, # Default values
                )

    with pytest.raises(
        ValueError,
        match=re.compile(
            pattern=(
                r"Can not provide multiple pseudos for the same element!.*"
                r"Problem detected for code 'pyscf'"
                ),
            flags=re.DOTALL,
            ),
        ):
            _ = PseudoSet.from_mixed_dir(
                pseudo_dir = psp_dir,
                codes      = {"pyscf"},
                extensions = None,
                )

    # Filter by include
    _ = PseudoSet.from_mixed_dir(
        pseudo_dir = psp_dir,
        codes      = {"pyscf"},
        include    = {"pyscf": "special"},
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
#end def test_from_mixed_dir_pyscf_collision


def test_from_mixed_dir_mismatches(tmp_path):
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
        match="Mismatch between provided extensions and codes!",
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
        match="Mismatch between provided include patterns and codes!",
        ):
        PseudoSet.from_mixed_dir(
            pseudo_dir = tmp_path,
            codes      = {"espresso"},
            include    = {
                "espresso": "*BFD*",
                "rmg": "*BFD*", # Not in the specified codes
                },
            )

    with pytest.raises(
        ValueError,
        match="Mismatch between provided exclude patterns and codes!",
        ):
        PseudoSet.from_mixed_dir(
            pseudo_dir = tmp_path,
            codes      = {"espresso"},
            exclude    = {
                "espresso": "*BFD*",
                "rmg": "*BFD*", # Not in the specified codes
                },
            )
#end def test_from_mixed_dir_mismatches


def test_from_mixed_dir_arg_for_all(tmp_path):
    pseudo_names = (
        "C.ccECP.xml",
        "H.ccECP.xml",
        "O.ccECP.xml",
        "C.ccECP.upf",
        "H.ccECP.upf",
        "O.ccECP.upf",
        "C.BFD.upf",
        "H.BFD.upf",
        "O.BFD.upf",
        "Fe.BFD.upf",
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
                "C": (psp_dir / "C.ccECP.xml").resolve(),
                "H": (psp_dir / "H.ccECP.xml").resolve(),
                "O": (psp_dir / "O.ccECP.xml").resolve(),
                },
            codes = {"qmcpack"},
            ),
        "espresso": PseudoSet(
            pseudos = {
                "C" : (psp_dir / "C.ccECP.upf").resolve(),
                "H" : (psp_dir / "H.ccECP.upf").resolve(),
                "O" : (psp_dir / "O.ccECP.upf").resolve(),
                },
            codes={"espresso"},
            ),
        "rmg": PseudoSet(
            pseudos = {
                "C": (psp_dir / "C.ccECP.xml").resolve(),
                "H": (psp_dir / "H.ccECP.xml").resolve(),
                "O": (psp_dir / "O.ccECP.xml").resolve(),
                },
            codes={"rmg"},
            ),
        }

    pseudoset = PseudoSet.from_mixed_dir(
        pseudo_dir = psp_dir,
        codes      = {"qmcpack", "rmg"},
        extensions = ".xml",
        )

    assert(pseudoset.keys() == {"qmcpack", "rmg"})

    qmcpack_calc = pseudoset["qmcpack"]
    qmcpack_ref = ref_pseudos["qmcpack"]
    assert(qmcpack_calc.pseudos     == qmcpack_ref.pseudos)
    assert(qmcpack_calc.codes       == qmcpack_ref.codes)
    assert(qmcpack_calc.pseudo_dirs == qmcpack_ref.pseudo_dirs)

    rmg_calc = pseudoset["rmg"]
    rmg_ref = ref_pseudos["rmg"]
    assert(rmg_calc.pseudos     == rmg_ref.pseudos)
    assert(rmg_calc.codes       == rmg_ref.codes)
    assert(rmg_calc.pseudo_dirs == rmg_ref.pseudo_dirs)

    pseudoset = PseudoSet.from_mixed_dir(
        pseudo_dir = psp_dir,
        codes      = {"qmcpack", "espresso"},
        include    = "*ccECP*",
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


    ref_Zeff_map = {
        "C": 4,
        "H": 1,
        "O": 6,
        }
    pseudoset = PseudoSet.from_mixed_dir(
        pseudo_dir    = psp_dir,
        codes         = {"qmcpack", "espresso"},
        include       = "*ccECP*",
        code_Zeff_map = ref_Zeff_map,
        )

    assert(set(pseudoset.keys()) == {"espresso", "qmcpack"})

    qmcpack_calc = pseudoset["qmcpack"]
    qmcpack_ref = ref_pseudos["qmcpack"]
    assert(qmcpack_calc.pseudos     == qmcpack_ref.pseudos)
    assert(qmcpack_calc.codes       == qmcpack_ref.codes)
    assert(qmcpack_calc.pseudo_dirs == qmcpack_ref.pseudo_dirs)
    assert(qmcpack_calc.Zeff_map    == ref_Zeff_map)

    espresso_calc = pseudoset["espresso"]
    espresso_ref = ref_pseudos["espresso"]
    assert(espresso_calc.pseudos     == espresso_ref.pseudos)
    assert(espresso_calc.codes       == espresso_ref.codes)
    assert(espresso_calc.pseudo_dirs == espresso_ref.pseudo_dirs)
    assert(espresso_calc.Zeff_map    == ref_Zeff_map)
#end def test_from_mixed_dir_arg_for_all


def test_priv_get_pseudos(tmp_path):
    qmcpack_dir, _, ref_qmcpack_pseudos = setup_psps(test_dir=tmp_path, code="qmcpack")
    qmcpack_pseudoset = PseudoSet.from_dir(pseudo_dir=qmcpack_dir, code ="detect")
    ref_pseudos = {
        ref_qmcpack_pseudos["C"].name: str(ref_qmcpack_pseudos["C"]),
        ref_qmcpack_pseudos["H"].name: str(ref_qmcpack_pseudos["H"]),
        }

    elements = ["C", "C", "H"]

    pseudos = qmcpack_pseudoset._get_pseudos(system=elements, code="qmcpack")

    assert(pseudos == ref_pseudos)

    system = generate_physical_system(
        elem = ["C", "H", "H", "H", "H"],
        pos = np.empty((5,3), dtype=float),
        C = 4,
        H = 1,
        )

    pseudos = qmcpack_pseudoset._get_pseudos(system=system, code="qmcpack")

    assert(pseudos == ref_pseudos)

    with pytest.raises(
        ValueError,
        match='Pseudopotential set is not available for code "espresso"'
        ):
        qmcpack_pseudoset._get_pseudos(system=elements, code="espresso")

    with pytest.raises(
        ValueError,
        match=r"Pseudopotential set does not contain the following species:\n.*Fe"
        ):
        qmcpack_pseudoset._get_pseudos(system=["C", "H", "Fe"], code="qmcpack")
#end def test_priv_get_pseudos


def test_get_Zeff():
    psp_dir = TEST_FILES["C.BFD.xml"].parent.resolve()
    ref_Zeff_default = {
        "C": 4,
        "H": 1,
        "O": 6,
        }

    qmcpack_pseudoset = PseudoSet.from_dir(
        pseudo_dir = psp_dir,
        code       = "qmcpack",
        extension  = None,
        )
    Zeff = qmcpack_pseudoset.get_Zeff(elem_labels=["C", "H", "O"])

    assert(Zeff == ref_Zeff_default)

    ref_Zeff_customized = {
        "C": 6,
        "H": 1,
        "O": 6,
        }

    qmcpack_pseudoset_custom = PseudoSet.from_dir(
        pseudo_dir = psp_dir,
        code       = "qmcpack",
        Zeff_map   = {"C": 6},
        extension  = None,
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
        pseudo_dir = psp_dir,
        code       = "qmcpack",
        extension  = None,
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


@pytest.mark.filterwarnings(r"ignore:.*ppset is deprecated.*")
@isolate_nexus_core
def test_legacy_ppset(tmp_path):
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

    ppset(
        label = 'shared_upf',
        pwscf = ['C.BFD.upf','H.BFD.upf','O.BFD.upf'],
        rmg   = ['C.BFD.upf','H.BFD.upf','O.BFD.upf'],
        )
    shared = PseudoSet.labeled_pseudosets['shared_upf']
    assert set(shared) == {'espresso','rmg'}
    assert shared['espresso'] is not shared['rmg']

    with pytest.raises(
        ValueError,
        match='label "bfd" is already registered',
        ):
        ppset(label='bfd',qmcpack=['C.BFD.xml','H.BFD.xml','O.BFD.xml'])

    with pytest.raises(
        ValueError,
        match='are not present in PseudoSet.pseudo_files',
        ):
        ppset(label='missing',qmcpack=['Ne.missing.xml'])

    with pytest.raises(ValueError, match='same elements'):
        ppset(
            label    = 'unequal_elements',
            espresso = ['C.BFD.upf','H.BFD.upf'],
            qmcpack  = ['C.BFD.xml','O.BFD.xml'],
            )

    with pytest.raises(ValueError,match='not compatible with that code'):
        ppset(label='incompatible',gamess=['C.BFD.xml','H.BFD.xml','O.BFD.xml'])
#end def test_legacy_ppset


@pytest.mark.filterwarnings(r"ignore:.*ppset is deprecated.*")
@isolate_nexus_core
def test_get_pseudos(tmp_path):
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
        label    = 'bfd',
        espresso = ["C.BFD.upf", "H.BFD.upf", "O.BFD.upf"],
        qmcpack  = ["C.BFD.xml", "H.BFD.xml", "O.BFD.xml"],
        gamess   = ["C.BFD.gms", "H.BFD.gms", "O.BFD.gms"],
        )

    system = generate_physical_system(
        elem = ['O','C','H'],
        pos  = np.empty((3,3),dtype=float),
        C    = 4,
        H    = 1,
        O    = 6,
        )
    remapped = PseudoSet.get_pseudos(
        pseudos = 'bfd',
        system = system,
        code = 'qmcpack',
        )
    assert list(remapped) == ['C.BFD.xml','H.BFD.xml','O.BFD.xml']
    assert remapped == {
        filename:PseudoSet.pseudo_files[filename] for filename in remapped
        }

    remapped = PseudoSet.get_pseudos(
        pseudos = ref_pseudos['bfd']['qmcpack'],
        system = system,
        code = 'qmcpack',
        )
    assert list(remapped) == ['C.BFD.xml','H.BFD.xml','O.BFD.xml']
    assert remapped == {
        filename:PseudoSet.pseudo_files[filename] for filename in remapped
        }

    explicit = ['C.BFD.xml', 'H.BFD.xml', 'O.BFD.xml']
    remapped = PseudoSet.get_pseudos(
        pseudos = explicit,
        system = system,
        code = 'qmcpack',
        )
    assert list(remapped) == explicit

    qmcpack_pseudoset = PseudoSet.labeled_pseudosets['bfd']['qmcpack']
    generated = {'qmcpack': qmcpack_pseudoset}
    with pytest.raises(
        TypeError,
        match='must contain only PseudoSet values',
        ):
        PseudoSet.get_pseudos(
            pseudos = {'qmcpack':qmcpack_pseudoset,'espresso':None},
            system = system,
            code = 'qmcpack',
            )

    with pytest.raises(
        ValueError,
        match='not available for code "gamess"',
        ):
        PseudoSet.get_pseudos(
            pseudos = generated,
            system = system,
            code = 'gamess',
            )

    with pytest.raises(
        ValueError,
        match='Either provide a list of pseudo names or use a `PhysicalSystem` object',
        ):
        PseudoSet.get_pseudos(
            pseudos = generated,
            system = None,
            code = 'qmcpack',
            )

    with pytest.raises(
        KeyError,
        match='label "unknown" is not registered',
        ):
        PseudoSet.get_pseudos(
            pseudos = 'unknown',
            system = system,
            code = 'qmcpack',
            )

    with pytest.raises(
        FileNotFoundError,
        match='are not present in PseudoSet.pseudo_files',
        ):
        PseudoSet.get_pseudos(
            pseudos = ['Ne.missing.xml'],
            system = system,
            code = 'qmcpack',
            )

    missing_species = generate_physical_system(
        elem = ['C','N'],
        pos  = np.empty((2,3),dtype=float),
        C    = 4,
        N    = 5,
        )
    with pytest.raises(
        ValueError,
        match=re.compile(
            pattern=r"Pseudopotential set does not contain the following species:.*\['N'\]",
            flags=re.DOTALL,
            )
        ):
        PseudoSet.get_pseudos(
            pseudos = generated,
            system = missing_species,
            code = 'qmcpack',
            )

    with pytest.raises(
        ValueError,
        match=re.compile(
            pattern=r"Pseudopotential set does not contain the following species:.*\['N'\]",
            flags=re.DOTALL,
            )
        ):
        PseudoSet.get_pseudos(
            pseudos = 'bfd',
            system = missing_species,
            code = 'qmcpack',
            )
#end def test_get_pseudos


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
)"""
    assert(repr(pseudoset) == ref_repr)
#end def test_pseudoset_repr


def test_from_mixed_dir_all_codes_partial_extensions(tmp_path):
    """Every explicitly selected code should be represented in the result."""
    (tmp_path / "C.upf").touch()

    pseudosets = PseudoSet.from_mixed_dir(
        pseudo_dir = tmp_path,
        codes      = {"qmcpack", "espresso"},
        extensions = {"qmcpack": {".xml"}},
        )

    assert(set(pseudosets) == {"qmcpack", "espresso"})
#end def test_from_mixed_dir_all_codes_partial_extensions


def test_normalize_code_map_keys():
    with pytest.raises(
        ValueError,
        match="Dictionary supplied two aliases for the same code, found duplicate:"
        ):
        PseudoSet._normalize_code_map_keys(
            {
                "pwscf": ("C.first.upf",),
                "espresso": ("C.second.upf",),
            }
        )
#end def test_normalize_code_map_keys


def test_from_mixed_dir_normalizes_code_alias_maps(tmp_path):
    """Aliases must be normalized consistently across all code-keyed inputs."""
    (tmp_path / "C.BFD.upf").touch()
    (tmp_path / "C.ONCV.upf").touch()

    pseudosets = PseudoSet.from_mixed_dir(
        pseudo_dir    = tmp_path,
        codes         = "pwscf",
        include       = {"pwscf": "*ONCV*"},
        code_Zeff_map = {"pwscf": {"C": 4}},
        )

    assert(set(pseudosets) == {"espresso"})
    assert(pseudosets["espresso"].pseudos == {
        "C": (tmp_path / "C.ONCV.upf").resolve()
        })
    assert(pseudosets["espresso"].Zeff_map == {"C": 4})
#end def test_from_mixed_dir_normalizes_code_alias_maps


def test_from_mixed_dir_returns_canonical_code_keys(tmp_path):
    """Result keys should use the same canonical lowercase code vocabulary."""
    (tmp_path / "C.xml").touch()

    pseudosets = PseudoSet.from_mixed_dir(
        pseudo_dir = tmp_path,
        codes      = {"QMCPACK"},
        extensions = {"QMCPACK": {".xml"}},
        )

    assert(set(pseudosets) == {"qmcpack"})
#end def test_from_mixed_dir_returns_canonical_code_keys


def test_from_mixed_dir_does_not_mutate_extensions(tmp_path):
    """Input mappings should remain unchanged."""
    extensions = {"qmcpack": {".xml"}}
    reference = deepcopy(extensions)

    PseudoSet.from_mixed_dir(
        pseudo_dir = tmp_path,
        extensions = extensions,
        )

    assert(extensions == reference)
#end def test_from_mixed_dir_does_not_mutate_extensions


def test_from_mixed_dir_accepts_immutable_extensions_mapping(tmp_path):
    """Any object satisfying the annotated ``Mapping`` contract should work."""
    extensions = MappingProxyType({"qmcpack": {".xml"}})

    pseudosets = PseudoSet.from_mixed_dir(
        pseudo_dir = tmp_path,
        extensions = extensions,
        )

    assert("qmcpack" in pseudosets)
#end def test_from_mixed_dir_accepts_immutable_extensions_mapping


def test_from_dir_ignores_potcar_directories(tmp_path):
    """A POTCAR must be a file, not merely an existing path."""
    element_dir = tmp_path / "C"
    element_dir.mkdir()
    (element_dir / "POTCAR").mkdir()
    with pytest.raises(NotADirectoryError, match="POTCARs can not be directories!"):
        PseudoSet.from_dir(
            pseudo_dir = tmp_path,
            code       = "vasp",
            )
#end def test_from_dir_ignores_potcar_directories


@isolate_nexus_core
def test_generate_pseudoset(tmp_path):
    pseudo_names = (
        "C.ccECP.gamess",
        "C.ccECP.nwchem",
        "C.ccECP.upf",
        "C.ccECP.xml",
        "H.ccECP.gamess",
        "H.ccECP.nwchem",
        "H.ccECP.upf",
        "H.ccECP.xml",
        "O.ccECP.gamess",
        "O.ccECP.nwchem",
        "O.ccECP.upf",
        "O.ccECP.xml",
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
        'espresso': PseudoSet(
            codes = {'espresso'},
            pseudos = {
                'C': (psp_dir / 'C.ccECP.upf').resolve(),
                'H': (psp_dir / 'H.ccECP.upf').resolve(),
                'O': (psp_dir / 'O.ccECP.upf').resolve(),
            },
            Zeff_map = {},
            ),
        'gamess': PseudoSet(
            codes = {'gamess'},
            pseudos = {
                'C': (psp_dir / 'C.ccECP.gamess').resolve(),
                'H': (psp_dir / 'H.ccECP.gamess').resolve(),
                'O': (psp_dir / 'O.ccECP.gamess').resolve(),
            },
            Zeff_map = {},
            ),
        'qmcpack': PseudoSet(
            codes = {'qmcpack'},
            pseudos = {
                'C': (psp_dir / 'C.ccECP.xml').resolve(),
                'H': (psp_dir / 'H.ccECP.xml').resolve(),
                'O': (psp_dir / 'O.ccECP.xml').resolve(),
            },
            Zeff_map = {},
            ),
        'pyscf': PseudoSet(
            codes = {'pyscf'},
            pseudos = {
                'C': (psp_dir / 'C.ccECP.nwchem').resolve(),
                'H': (psp_dir / 'H.ccECP.nwchem').resolve(),
                'O': (psp_dir / 'O.ccECP.nwchem').resolve(),
            },
            Zeff_map = {},
            ),
        'rmg': PseudoSet(
            codes = {'rmg'},
            pseudos = {
                'C': (psp_dir / 'C.ccECP.xml').resolve(),
                'H': (psp_dir / 'H.ccECP.xml').resolve(),
                'O': (psp_dir / 'O.ccECP.xml').resolve(),
            },
            Zeff_map = {},
            )
        }

    psps = generate_pseudoset(
        pseudo_dir=psp_dir,
        code={"qmcpack", "espresso", "rmg", "pyscf", "gamess"},
        extension={"rmg": ".xml"},
        )

    assert(psps.keys() == ref_pseudos.keys())
    for code, ref_psp_set in ref_pseudos.items():
        gen_psp_set = psps[code]
        assert(gen_psp_set.pseudos == ref_psp_set.pseudos)
        assert(gen_psp_set.codes == ref_psp_set.codes)
        assert(gen_psp_set.Zeff_map == ref_psp_set.Zeff_map)
        assert(gen_psp_set.pseudo_dirs == ref_psp_set.pseudo_dirs)

    # Test ppset-like interface
    psps = generate_pseudoset(
        pseudo_dir=psp_dir,
        qmcpack  = ["C.ccECP.xml",    "H.ccECP.xml",    "O.ccECP.xml"],
        espresso = ["C.ccECP.upf",    "H.ccECP.upf",    "O.ccECP.upf"],
        rmg      = ["C.ccECP.xml",    "H.ccECP.xml",    "O.ccECP.xml"],
        pyscf    = ["C.ccECP.nwchem", "H.ccECP.nwchem", "O.ccECP.nwchem"],
        gamess   = ["C.ccECP.gamess", "H.ccECP.gamess", "O.ccECP.gamess"],
        )

    assert(psps.keys() == ref_pseudos.keys())
    for code, ref_psp_set in ref_pseudos.items():
        gen_psp_set = psps[code]
        assert(gen_psp_set.pseudos == ref_psp_set.pseudos)
        assert(gen_psp_set.Zeff_map == ref_psp_set.Zeff_map)
        assert(gen_psp_set.pseudo_dirs == ref_psp_set.pseudo_dirs)
        # Switch to superset because codes are detected dynamically, not set statically
        assert(gen_psp_set.codes >= ref_psp_set.codes)

    # Test ppset-like interface, but with directories

    qmcpack_dir,  _, ref_qmcpack_pseudos  = setup_psps(test_dir=tmp_path, code="qmcpack")
    espresso_dir, _, ref_espresso_pseudos = setup_psps(test_dir=tmp_path, code="espresso")
    gamess_dir,   _, ref_gamess_pseudos   = setup_psps(test_dir=tmp_path, code="gamess")
    rmg_dir,      _, ref_rmg_pseudos      = setup_psps(test_dir=tmp_path, code="rmg")
    pyscf_dir,    _, ref_pyscf_pseudos    = setup_psps(test_dir=tmp_path, code="pyscf")

    ref_pseudos = {
        "qmcpack" : PseudoSet(ref_qmcpack_pseudos),
        "espresso": PseudoSet(ref_espresso_pseudos),
        "gamess"  : PseudoSet(ref_gamess_pseudos),
        "rmg"     : PseudoSet(ref_rmg_pseudos),
        "pyscf"   : PseudoSet(ref_pyscf_pseudos),
        }

    psps = generate_pseudoset(
        qmcpack  = qmcpack_dir,
        espresso = espresso_dir,
        rmg      = rmg_dir,
        pyscf    = pyscf_dir,
        gamess   = gamess_dir,
        )

    assert(psps.keys() == ref_pseudos.keys())
    for code, ref_psp_set in ref_pseudos.items():
        gen_psp_set = psps[code]
        assert(gen_psp_set.pseudos == ref_psp_set.pseudos)
        assert(gen_psp_set.codes == ref_psp_set.codes)
        assert(gen_psp_set.Zeff_map == ref_psp_set.Zeff_map)
        assert(gen_psp_set.pseudo_dirs == ref_psp_set.pseudo_dirs)

    # Test pseudo dirs relative to a path
    psps = generate_pseudoset(
        pseudo_dir = tmp_path,
        qmcpack  = qmcpack_dir.name,
        espresso = espresso_dir.name,
        rmg      = rmg_dir.name,
        pyscf    = pyscf_dir.name,
        gamess   = gamess_dir.name,
        )

    assert(psps.keys() == ref_pseudos.keys())
    for code, ref_psp_set in ref_pseudos.items():
        gen_psp_set = psps[code]
        assert(gen_psp_set.pseudos == ref_psp_set.pseudos)
        assert(gen_psp_set.codes == ref_psp_set.codes)
        assert(gen_psp_set.Zeff_map == ref_psp_set.Zeff_map)
        assert(gen_psp_set.pseudo_dirs == ref_psp_set.pseudo_dirs)

    with pytest.raises(
        ValueError,
        match="Must supply `pseudo_dir` and/or `codes_psps`!",
        ):
        generate_pseudoset()

    with pytest.raises(
        FileNotFoundError,
        match="`pseudo_dir` must exist!" ,
        ):
        generate_pseudoset(pseudo_dir="/path/to/nowhere")

    with pytest.raises(
        NotADirectoryError,
        match="`pseudo_dir` must be a directory!" ,
        ):
        generate_pseudoset(pseudo_dir=ref_qmcpack_pseudos["C"])

    with pytest.raises(
        ValueError,
        match="When supplying a direct map of codes to pseudos you cannot pass `code`!",
        ):
        generate_pseudoset(qmcpack=[ref_qmcpack_pseudos["C"]], code="qmcpack")

    with pytest.raises(
        TypeError,
        match="Must supply a directory or collection of file paths for direct map!",
        ):
        generate_pseudoset(qmcpack=1234)

    with pytest.raises(
        FileNotFoundError,
        match="The path for code qmcpack does not exist!",
        ):
        generate_pseudoset(qmcpack="/path/to/nowhere")

    with pytest.raises(
        RuntimeError,
        match="Error when processing pseudo directory for code 'qmcpack'",
        ):
        generate_pseudoset(qmcpack=psp_dir)

    with pytest.raises(
        NotADirectoryError,
        match="If you are providing a single path it must be to a directory!",
        ):
        generate_pseudoset(qmcpack=ref_qmcpack_pseudos["C"])
#end def test_generate_pseudoset
