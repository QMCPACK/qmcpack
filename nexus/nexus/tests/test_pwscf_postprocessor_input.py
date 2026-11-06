import pytest

from ..developer import obj, to_obj
from ..pwscf_postprocessors import (
    BandsInput,
    CpppInput,
    DosInput,
    PPInput,
    ProjwfcInput,
    PwexportInput,
    generate_bands_input,
    generate_cppp_input,
    generate_dos_input,
    generate_pp_input,
    generate_projwfc_input,
    generate_pwexport_input,
)
from ..testing import object_eq
from . import NexusTestOrder

pytestmark = pytest.mark.order(NexusTestOrder.PWSCF_POSTPROCESSOR_INPUT)

projwfc_in = '''&projwfc
  prefix = 'pwscf'
  outdir = 'pwscf_output'
/
'''

def test_empty_init():
    ppi = PPInput()
    ppi = generate_pp_input()

    ppi = DosInput()
    ppi = generate_dos_input()

    ppi = BandsInput()
    ppi = generate_bands_input()

    ppi = ProjwfcInput()
    ppi = generate_projwfc_input()

    ppi = CpppInput()
    ppi = generate_cppp_input()

    ppi = PwexportInput()
    ppi = generate_pwexport_input()
#end def test_empty_init



def test_read(tmp_path):
    infile_path = tmp_path / 'projwfc.in'
    infile_path.write_text(projwfc_in)

    pi = ProjwfcInput(infile_path)

    pi_ref = obj(
        projwfc = obj(
            prefix = 'pwscf',
            outdir = 'pwscf_output',
            ),
        )

    assert(object_eq(to_obj(pi),pi_ref))
#end def test_read



def test_write(tmp_path):
    infile_path = tmp_path / 'projwfc.in'
    infile_path.write_text(projwfc_in)

    write_path = tmp_path / 'projwfc_write.in'
    pi_write = ProjwfcInput(infile_path)

    pi_write.write(write_path)

    pi_read = ProjwfcInput(write_path)

    pi_ref = obj(
        projwfc = obj(
            prefix = 'pwscf',
            outdir = 'pwscf_output',
            ),
        )

    assert(object_eq(to_obj(pi_read),pi_ref))
#end def test_write



def test_generate():
    pi = generate_projwfc_input(
        prefix = 'pwscf',
        outdir = 'pwscf_output',
        )

    pi_ref = obj(
        projwfc = obj(
            prefix = 'pwscf',
            outdir = 'pwscf_output',
            ),
        )

    assert(object_eq(to_obj(pi),pi_ref))
#end def test_generate
