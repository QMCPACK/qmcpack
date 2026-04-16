try:
    import pytest
    from . import NexusTestOrder
    pytestmark = pytest.mark.order(NexusTestOrder.USER_EXAMPLES)
except ImportError:
    pass

from pathlib import Path
import os
import sys
import shutil
from subprocess import Popen, PIPE


nexus_root = Path(__file__).resolve().parent.parent.parent # qmcpack/nexus
example_root  = nexus_root / "nexus/examples"
test_root     = nexus_root / "nexus/tests"
reference_dir = test_root / "reference/user_examples"
output_root   = test_root / "test_user_examples_output"

qmcpack_pseudos  = example_root / "qmcpack/pseudopotentials"
espresso_pseudos = example_root / "quantum_espresso/pseudopotentials"


def copy_pseudos(code: str):

    if code == "qmcpack":
        output_path = output_root / "qmcpack/pseudopotentials"
        shutil.copytree(qmcpack_pseudos, output_path, dirs_exist_ok=True)
    elif code == "quantum_espresso":
        output_path = output_root / "quantum_espresso/pseudopotentials"
        shutil.copytree(espresso_pseudos, output_path, dirs_exist_ok=True)
    else:
        raise ValueError(f"Invalid code for pseudopotential identification: {code}")


def copy_example_files(example_dir: str):

    example_path = example_root / example_dir
    output_path = output_root / example_dir
    test_path = shutil.copytree(example_path, output_path, dirs_exist_ok=True)

    return test_path


def run_example_script(script: str, test_path: Path):

    old_cwd = Path.cwd()
    os.chdir(test_path)

    script_command = f"PYTHONPATH={nexus_root} {sys.executable} ./{script} --generate_only --sleep=0.01"
    process = Popen(script_command, shell=True, stdout=PIPE, stderr=PIPE, close_fds=True)
    out, err = process.communicate()
    returncode = process.returncode
    os.chdir(old_cwd) # Reset us back to the old cwd
    if returncode != 0:
        msg = (
            "Executed system command failed.\n"
            "\n"
            "Command:\n"
            "========\n"
           f"{script_command}\n"
            "\n"
            "stdout:\n"
            "=======\n"
           f"{out}\n"
            "\n"
            "stderr:\n"
            "=======\n"
           f"{err}\n"
            "\n"
            "Return code:\n"
            "============\n"
           f"{returncode}\n"
        )
        # Return test failure so it can be dealt with in the respective test.
        return False, msg
    else:
        return True, "Success!"


def check_generated_files(
    code: str,
    example_path: str,
    filepath: str,
):

    from ..testing import object_diff
    from nexus.pwscf_input import PwscfInput
    from nexus.qmcpack_converters import Pw2qmcpackInput
    from nexus.gamess_input import GamessInput
    from nexus.qmcpack_input import QmcpackInput

    input_classes = dict(
        pwscf      = PwscfInput,
        pw2qmcpack = Pw2qmcpackInput,
        gamess     = GamessInput,
        qmcpack    = QmcpackInput,
    )

    ref_filepath = reference_dir / example_path / filepath
    gen_filepath = output_root / example_path / filepath

    if not ref_filepath.exists():
        raise FileNotFoundError(
            "Reference file is missing\n"
            f"File should be located at: {ref_filepath!s}"
        )
    elif not gen_filepath.exists():
        raise FileNotFoundError(
            f"Input file was not generated: {gen_filepath!s}"
        )

    input_class = input_classes[code]
    ref_input = input_class(str(ref_filepath))
    gen_input = input_class(str(gen_filepath))
    diff, dgen, dref = object_diff(gen_input, ref_input, full=True)
    if diff:
        # assume failure
        failed = True
        # if difference due to periodically equivalent points
        # then it is not a failure
        check_pbc = False
        if len(dgen)==1 and len(dref)==1:
            kgen = list(dgen.keys())[0].rsplit('/',1)[1]
            kref = list(dref.keys())[0].rsplit('/',1)[1]
            check_pbc |= code=='qmcpack' and kgen==kref=='position'
            check_pbc |= code=='pwscf'   and kgen==kref=='positions'
        #end if
        if check_pbc:
            try:
                # extract Structure objects from SimulationInput objects
                rs = ref_input.return_structure()
                gs = gen_input.return_structure()
                # compare minimum image distances of all atomic coordinates
                d = rs.min_image_distances(gs.pos,pairs=False)
                # allow for small deviation due to precision of ascii floats in the text input files 
                if d.min()<1e-6:
                    failed = False
                #end if
            except:
                None
            #end try
        #end if
        if failed:
            # report on failures
            from nexus.generic import obj
            dgen = obj(dgen)
            dref = obj(dref)
            msg  = 'reference and generated input files differ\n'
            msg += 'reference file: '+os.path.realpath(filepath)+'\n'
            msg += 'reference file difference\n'
            msg += 40*'='+'\n'
            msg += str(dref)
            msg += 'generated file difference\n'
            msg += 40*'='+'\n'
            msg += str(dgen)
            return False, msg
        #end if
    #end if
    return True, "Success!"


# Move pseudos
copy_pseudos("qmcpack")
copy_pseudos("quantum_espresso")


def test_pwscf_relax_Ge_T():

    test_data = dict(
        path = 'quantum_espresso/relax_Ge_T_vs_kpoints', 
        scripts = [
            'relax_vs_kpoints_example.py',
        ],
        files = [
            ('pwscf', 'input', 'runs/relax/kgrid_111/relax.in'),
            ('pwscf', 'input', 'runs/relax/kgrid_222/relax.in'),
            ('pwscf', 'input', 'runs/relax/kgrid_444/relax.in'),
            ('pwscf', 'input', 'runs/relax/kgrid_666/relax.in'),
        ],
    )

    test_path = copy_example_files(test_data["path"])

    for script in test_data["scripts"]:
        success, message = run_example_script(script, test_path)
        assert(success), message

    for code, filetype, filepath in test_data["files"]:
        success, message = check_generated_files(
            code,
            test_data["path"],
            filepath,
        )
        assert(success), message


def test_gamess_H2O():

    test_data = dict(
        path = 'gamess/H2O',
        scripts = [
            'h2o_pp_hf.py',
            'h2o_pp_cisd.py',
            'h2o_pp_casscf.py',
        ],
        files = [
            ('gamess', 'input', 'runs/pp_hf/rhf.inp'),
            ('gamess', 'input', 'runs/pp_cisd/rhf.inp'),
            ('gamess', 'input', 'runs/pp_cisd/cisd.inp'),
            ('gamess', 'input', 'runs/pp_casscf/rhf.inp'),
            ('gamess', 'input', 'runs/pp_casscf/cas.inp'),
        ],
    )

    test_path = copy_example_files(test_data["path"])

    for script in test_data["scripts"]:
        success, message = run_example_script(script, test_path)
        assert(success), message

    for code, filetype, filepath in test_data["files"]:
        success, message = check_generated_files(
            code,
            test_data["path"],
            filepath,
        )
        assert(success), message


def test_qmcpack_H2O():

    test_data = dict(
        path = 'qmcpack/rsqmc_misc/H2O',
        scripts = [
            'H2O.py',
        ],
        files = [
            ('pwscf',      'input', 'runs/scf.in'),
            ('pw2qmcpack', 'input', 'runs/p2q.in'),
            ('qmcpack',    'input', 'runs/opt.in.xml'),
            ('qmcpack',    'input', 'runs/dmc.in.xml'),
        ],
    )

    test_path = copy_example_files(test_data["path"])

    for script in test_data["scripts"]:
        success, message = run_example_script(script, test_path)
        assert(success), message

    for code, filetype, filepath in test_data["files"]:
        success, message = check_generated_files(
            code,
            test_data["path"],
            filepath,
        )
        assert(success), message


def test_qmcpack_LiH():

    test_data = dict(
        path = 'qmcpack/rsqmc_misc/LiH',
        scripts = [
            'LiH.py',
        ],
        files = [
            ('pwscf',      'input', 'runs/scf.in'),
            ('pwscf',      'input', 'runs/nscf.in'),
            ('pw2qmcpack', 'input', 'runs/p2q.in'),
            ('qmcpack',    'input', 'runs/opt.in.xml'),
            ('qmcpack',    'input', 'runs/dmc.in.xml'),
        ],
    )

    test_path = copy_example_files(test_data["path"])

    for script in test_data["scripts"]:
        success, message = run_example_script(script, test_path)
        assert(success), message

    for code, filetype, filepath in test_data["files"]:
        success, message = check_generated_files(
            code,
            test_data["path"],
            filepath,
        )
        assert(success), message


def test_qmcpack_c20():

    test_data = dict(
        path = 'qmcpack/rsqmc_misc/c20',
        scripts = [
            'c20.py',
        ],
        files = [
            ('pwscf',      'input', 'runs/c20/scf/scf.in'),
            ('pw2qmcpack', 'input', 'runs/c20/nscf/p2q.in'),
            ('qmcpack',    'input', 'runs/c20/opt/opt.in.xml'),
            ('qmcpack',    'input', 'runs/c20/qmc/qmc.in.xml'),
        ],
    )

    test_path = copy_example_files(test_data["path"])

    for script in test_data["scripts"]:
        success, message = run_example_script(script, test_path)
        assert(success), message

    for code, filetype, filepath in test_data["files"]:
        success, message = check_generated_files(
            code,
            test_data["path"],
            filepath,
        )
        assert(success), message


def test_qmcpack_diamond():

    test_data = dict(
        path = 'qmcpack/rsqmc_misc/diamond',
        scripts = [
            'diamond.py',
            'diamond_vacancy.py',
        ],
        files = [
            ('pwscf',      'input', 'runs/diamond/scf/scf.in'),
            ('pw2qmcpack', 'input', 'runs/diamond/scf/conv.in'),
            ('qmcpack',    'input', 'runs/diamond/vmc/vmc.in.xml'),
            ('pwscf',      'input', 'runs/diamond_vacancy/relax/relax.in'),
            ('pwscf',      'input', 'runs/diamond_vacancy/scf/scf.in'),
        ],
    )

    test_path = copy_example_files(test_data["path"])

    for script in test_data["scripts"]:
        success, message = run_example_script(script, test_path)
        assert(success), message

    for code, filetype, filepath in test_data["files"]:
        success, message = check_generated_files(
            code,
            test_data["path"],
            filepath,
        )
        assert(success), message


def test_qmcpack_graphene():

    test_data = dict(
        path = 'qmcpack/rsqmc_misc/graphene',
        scripts = [
            'graphene.py',
        ],
        files = [
            ('pwscf',      'input', 'runs/graphene/scf/scf.in'),
            ('pwscf',      'input', 'runs/graphene/nscf/nscf.in'),
            ('pw2qmcpack', 'input', 'runs/graphene/nscf/p2q.in'),
            ('pwscf',      'input', 'runs/graphene/nscf_opt/nscf.in'),
            ('pw2qmcpack', 'input', 'runs/graphene/nscf_opt/p2q.in'),
            ('qmcpack',    'input', 'runs/graphene/opt/opt.in.xml'),
            ('qmcpack',    'input', 'runs/graphene/qmc/qmc.in.xml'),
        ],
    )

    test_path = copy_example_files(test_data["path"])

    for script in test_data["scripts"]:
        success, message = run_example_script(script, test_path)
        assert(success), message

    for code, filetype, filepath in test_data["files"]:
        success, message = check_generated_files(
            code,
            test_data["path"],
            filepath,
        )
        assert(success), message


def test_qmcpack_oxygen_dimer():

    test_data = dict(
        path = 'qmcpack/rsqmc_misc/oxygen_dimer',
        scripts = [
            'oxygen_dimer.py',
        ],
        files = [
            ('pwscf',      'input', 'scale_1.0/scf.in'),
            ('pw2qmcpack', 'input', 'scale_1.0/p2q.in'),
            ('qmcpack',    'input', 'scale_1.0/opt.in.xml'),
            ('qmcpack',    'input', 'scale_1.0/qmc.in.xml'),
        ],
    )

    test_path = copy_example_files(test_data["path"])

    for script in test_data["scripts"]:
        success, message = run_example_script(script, test_path)
        assert(success), message

    for code, filetype, filepath in test_data["files"]:
        success, message = check_generated_files(
            code,
            test_data["path"],
            filepath,
        )
        assert(success), message
