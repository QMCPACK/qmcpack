
import testing
from testing import value_eq,object_eq,text_eq


projwfc_in = '''&projwfc
  prefix = 'pwscf'
  outdir = 'pwscf_output'
/
'''

pw2wan90_in = '''&inputpp
  prefix = 'pwscf'
  outdir = 'pwscf_output'
  seedname = 'wannier'
/
'''


def test_import():
    from pwscf_postprocessors import PPInput,generate_pp_input
    from pwscf_postprocessors import DosInput,generate_dos_input
    from pwscf_postprocessors import BandsInput,generate_bands_input
    from pwscf_postprocessors import ProjwfcInput,generate_projwfc_input
    from pwscf_postprocessors import CpppInput,generate_cppp_input
    from pwscf_postprocessors import PwexportInput,generate_pwexport_input
    from pwscf_postprocessors import Pw2wan90Input,generate_pw2wan90_input
#end def test_import



def test_empty_init():
    from pwscf_postprocessors import PPInput,generate_pp_input
    from pwscf_postprocessors import DosInput,generate_dos_input
    from pwscf_postprocessors import BandsInput,generate_bands_input
    from pwscf_postprocessors import ProjwfcInput,generate_projwfc_input
    from pwscf_postprocessors import CpppInput,generate_cppp_input
    from pwscf_postprocessors import PwexportInput,generate_pwexport_input
    from pwscf_postprocessors import Pw2wan90Input,generate_pw2wan90_input
    
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

    ppi = Pw2wan90Input()
    ppi = generate_pw2wan90_input() 

#end def test_empty_init



def test_read():
    import os
    from developer import obj
    from pwscf_postprocessors import ProjwfcInput

    tpath = testing.setup_unit_test_output_directory('pwscf_postprocessor_input','test_read')

    infile_path = os.path.join(tpath,'projwfc.in')
    open(infile_path,'w').write(projwfc_in)

    pi = ProjwfcInput(infile_path)
    
    pi_ref = obj(
        projwfc = obj(
            prefix = 'pwscf',
            outdir = 'pwscf_output',
            ),
        )

    assert(object_eq(pi.to_obj(),pi_ref))
#end def test_read



def test_write():
    import os
    from developer import obj
    from pwscf_postprocessors import ProjwfcInput

    tpath = testing.setup_unit_test_output_directory('pwscf_postprocessor_input','test_write')

    infile_path = os.path.join(tpath,'projwfc.in')
    open(infile_path,'w').write(projwfc_in)

    write_path = os.path.join(tpath,'projwfc_write.in')
    pi_write = ProjwfcInput(infile_path)
    
    pi_write.write(write_path)

    pi_read = ProjwfcInput(write_path)

    pi_ref = obj(
        projwfc = obj(
            prefix = 'pwscf',
            outdir = 'pwscf_output',
            ),
        )

    assert(object_eq(pi_read.to_obj(),pi_ref))
#end def test_write



def test_generate():
    from developer import obj
    from pwscf_postprocessors import generate_projwfc_input

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

    assert(object_eq(pi.to_obj(),pi_ref))
#end def test_generate

def test_pw2wan90_generate():
    from generic import obj
    from pwscf_postprocessors import generate_pw2wan90_input

    pi = generate_pw2wan90_input(
        prefix = 'pwscf',
        outdir = 'pwscf_output',
        seedname = 'wannier',
        )

    pi_ref = obj(
        inputpp = obj(
            prefix = 'pwscf',
            outdir = 'pwscf_output',
            seedname = 'wannier',
            ),
        )

    assert(object_eq(pi.to_obj(),pi_ref))   


def test_read_pw2wan90_input():
    import os
    from generic import obj
    from pwscf_postprocessors import Pw2wan90Input

    tpath = testing.setup_unit_test_output_directory('pwscf_postprocessor_input','test_read')

    infile_path = os.path.join(tpath,'pw2wan90.in')
    open(infile_path,'w').write(pw2wan90_in)

    pi = Pw2wan90Input(infile_path)
    
    pi_ref = obj(
        inputpp = obj(
            prefix = 'pwscf',
            outdir = 'pwscf_output',
            seedname = 'wannier',
            ),
        )

    assert(object_eq(pi.to_obj(),pi_ref))
#end def test_read



def test_write_pw2wan90_input():
    import os
    from generic import obj
    from pwscf_postprocessors import Pw2wan90Input

    tpath = testing.setup_unit_test_output_directory('pwscf_postprocessor_input','test_write')

    infile_path = os.path.join(tpath,'pw2wan90.in')
    open(infile_path,'w').write(pw2wan90_in)

    write_path = os.path.join(tpath,'pw2wan90_write.in')
    pi_write = Pw2wan90Input(infile_path)
    
    pi_write.write(write_path)

    pi_read = Pw2wan90Input(write_path)

    pi_ref = obj(
        inputpp = obj(
            prefix = 'pwscf',
            outdir = 'pwscf_output',
            seedname = 'wannier',
            ),
        )

    assert(object_eq(pi_read.to_obj(),pi_ref))
#end def test_write

print(test_read_pw2wan90_input())
print(test_write_pw2wan90_input())
