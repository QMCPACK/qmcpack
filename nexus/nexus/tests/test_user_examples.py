import os
import shutil
import sys
from pathlib import Path
from shutil import ignore_patterns
from subprocess import PIPE, Popen

import pytest

from ..developer import obj
from ..generic import generic_settings
from ..testing import object_diff, text_diff
from ..pwscf_input import PwscfInput
from ..qmcpack_converters import Pw2qmcpackInput
from ..gamess_input import GamessInput
from ..qmcpack_input import QmcpackInput
from ..rmg_input import RmgInput
from ..pyscf_input import PyscfInput
from . import TEST_DIR, NexusTestOrder

pytestmark = pytest.mark.order(NexusTestOrder.USER_EXAMPLES)
generic_settings.raise_error = True

nexus_root = TEST_DIR.parent.parent # qmcpack/nexus
example_root  = nexus_root / "nexus/examples"
test_root     = nexus_root / "nexus/tests"
reference_dir = test_root / "reference/user_examples"

qmcpack_pseudos  = example_root / "qmcpack/pseudopotentials"
espresso_pseudos = example_root / "quantum_espresso/pseudopotentials"


def copy_pseudos(code: str, tmp_dir: Path):
    """Copy pseudopotential files over to the test directories."""
    if code == "qmcpack":
        output_path = tmp_dir / "qmcpack/pseudopotentials"
        shutil.copytree(qmcpack_pseudos, output_path, dirs_exist_ok=True)
    elif code == "quantum_espresso":
        output_path = tmp_dir / "quantum_espresso/pseudopotentials"
        shutil.copytree(espresso_pseudos, output_path, dirs_exist_ok=True)
    else:
        raise ValueError(f"Invalid code for pseudopotential identification: {code}")
#end def copy_pseudos


def copy_example_files(
    example_dir: str,
    tmp_dir: Path,
    *,
    template_file: Path | None = None
    ):

    example_path = example_root / example_dir
    output_path = tmp_dir / example_dir

    test_path = shutil.copytree(
        src           = example_path,
        dst           = output_path,
        dirs_exist_ok = True,
        ignore        = ignore_patterns("*.py", "*.rst"),
        )

    if template_file is not None:
        shutil.copy(template_file, test_path)

    return test_path
#end def copy_example_files


def run_example_script(script: Path, test_path: Path) -> tuple[bool, str]:

    script = script.resolve() # Absolute path helps in debugging
    old_cwd = Path.cwd()
    os.chdir(test_path)

    script_command = f"PYTHONPATH={nexus_root} {sys.executable} {script} --generate_only --sleep=0.01"
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
           f"{str(out, encoding='utf-8')}\n"
            "\n"
            "stderr:\n"
            "=======\n"
           f"{str(err, encoding='utf-8')}\n"
            "\n"
            "Return code:\n"
            "============\n"
           f"{returncode}\n"
            )
        # Return test failure so it can be dealt with in the respective test.
        return False, msg
    else:
        return True, "Success!"
#end def run_example_script


def check_generated_files(
    code: str,
    tmp_dir: Path,
    example_path: str,
    filepath: str,
    ) -> tuple[bool, str]:
    input_classes = dict(
        pwscf      = PwscfInput,
        pw2qmcpack = Pw2qmcpackInput,
        gamess     = GamessInput,
        qmcpack    = QmcpackInput,
        rmg        = RmgInput,
        pyscf      = PyscfInput,
        )

    ref_filepath = reference_dir / example_path / filepath
    gen_filepath = tmp_dir / example_path / filepath

    if not ref_filepath.exists():
        msg = (
            "Reference file is missing\n"
            f"File should be located at: {ref_filepath!s}"
            )
        raise FileNotFoundError(msg)
    elif not gen_filepath.exists():
        msg = f"Input file was not generated: {gen_filepath!s}"
        raise FileNotFoundError(msg)

    if code == "pyscf":
        # PyscfInput does template string things, which we can't really
        # evaluate, so instead we just do a direct text comparison of the
        # generated output from the test.
        ref_input = ref_filepath.read_text()
        gen_input = gen_filepath.read_text()
        diff, dgen, dref = text_diff(ref_input, gen_input, full=True)
    else:
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
                pass
            #end try
        #end if
        if failed:
            # report on failures
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
#end def check_generated_files


def test_pwscf_relax_Ge_T(tmp_path):
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

    test_path = copy_example_files(test_data["path"], tmp_path)
    copy_pseudos("quantum_espresso", tmp_path)

    for script in test_data["scripts"]:
        script_path = example_root / test_data["path"] / script
        success, message = run_example_script(script_path, test_path)
        assert(success), message

    for code, _filetype, filepath in test_data["files"]:
        success, message = check_generated_files(
            code,
            tmp_path,
            test_data["path"],
            filepath,
            )
        assert(success), message
#end def test_pwscf_relax_Ge_T


def test_gamess_H2O(tmp_path):
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

    test_path = copy_example_files(test_data["path"], tmp_path)

    for script in test_data["scripts"]:
        script_path = example_root / test_data["path"] / script
        success, message = run_example_script(script_path, test_path)
        assert(success), message

    for code, _filetype, filepath in test_data["files"]:
        success, message = check_generated_files(
            code,
            tmp_path,
            test_data["path"],
            filepath,
            )
        assert(success), message
#end def test_gamess_H2O


def test_qmcpack_H2O(tmp_path):
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

    test_path = copy_example_files(test_data["path"], tmp_path)
    copy_pseudos("quantum_espresso", tmp_path)
    copy_pseudos("qmcpack", tmp_path)


    for script in test_data["scripts"]:
        script_path = example_root / test_data["path"] / script
        success, message = run_example_script(script_path, test_path)
        assert(success), message

    for code, _filetype, filepath in test_data["files"]:
        success, message = check_generated_files(
            code,
            tmp_path,
            test_data["path"],
            filepath,
            )
        assert(success), message
#end def test_qmcpack_H2O


def test_qmcpack_LiH(tmp_path):
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

    test_path = copy_example_files(test_data["path"], tmp_path)
    copy_pseudos("quantum_espresso", tmp_path)
    copy_pseudos("qmcpack", tmp_path)

    for script in test_data["scripts"]:
        script_path = example_root / test_data["path"] / script
        success, message = run_example_script(script_path, test_path)
        assert(success), message

    for code, _filetype, filepath in test_data["files"]:
        success, message = check_generated_files(
            code,
            tmp_path,
            test_data["path"],
            filepath,
            )
        assert(success), message
#end def test_qmcpack_LiH


def test_qmcpack_c20(tmp_path):
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

    test_path = copy_example_files(test_data["path"], tmp_path)
    copy_pseudos("quantum_espresso", tmp_path)
    copy_pseudos("qmcpack", tmp_path)

    for script in test_data["scripts"]:
        script_path = example_root / test_data["path"] / script
        success, message = run_example_script(script_path, test_path)
        assert(success), message

    for code, _filetype, filepath in test_data["files"]:
        success, message = check_generated_files(
            code,
            tmp_path,
            test_data["path"],
            filepath,
            )
        assert(success), message
#end def test_qmcpack_c20


def test_qmcpack_diamond(tmp_path):
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

    test_path = copy_example_files(test_data["path"], tmp_path)
    copy_pseudos("quantum_espresso", tmp_path)
    copy_pseudos("qmcpack", tmp_path)

    for script in test_data["scripts"]:
        script_path = example_root / test_data["path"] / script
        success, message = run_example_script(script_path, test_path)
        assert(success), message

    for code, _filetype, filepath in test_data["files"]:
        success, message = check_generated_files(
            code,
            tmp_path,
            test_data["path"],
            filepath,
            )
        assert(success), message
#end def test_qmcpack_diamond


@pytest.mark.parametrize(
    "example_dir,script",
    [
        (
            '05_diamond_dft_dmc_restart',
            'diamond_lda_dmc_restart_same_dir.py',
            ),
        (
            '05_diamond_dft_dmc_restart',
            'diamond_lda_dmc_restart_separate_dirs.py',
            ),
        (
            '06_diamond_dft_dmc_twistavg_restart',
            'diamond_lda_dmc_twistavg_restart_same_dir.py',
            ),
        (
            '06_diamond_dft_dmc_twistavg_restart',
            'diamond_lda_dmc_twistavg_restart_separate_dirs.py',
            ),
        ],
    )
def test_qmcpack_restart_examples(tmp_path,example_dir,script):
    path = 'qmcpack/rsqmc_quantum_espresso/'+example_dir
    test_path = copy_example_files(path,tmp_path)
    copy_pseudos('qmcpack',tmp_path)

    script_path = example_root / path / script
    success,message = run_example_script(script_path,test_path)
    assert(success),message
#end def test_qmcpack_restart_examples


def test_qmcpack_graphene(tmp_path):
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

    test_path = copy_example_files(test_data["path"], tmp_path)
    copy_pseudos("quantum_espresso", tmp_path)
    copy_pseudos("qmcpack", tmp_path)

    for script in test_data["scripts"]:
        script_path = example_root / test_data["path"] / script
        success, message = run_example_script(script_path, test_path)
        assert(success), message

    for code, _filetype, filepath in test_data["files"]:
        success, message = check_generated_files(
            code,
            tmp_path,
            test_data["path"],
            filepath,
            )
        assert(success), message
#end def test_qmcpack_graphene


def test_qmcpack_oxygen_dimer(tmp_path):
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

    test_path = copy_example_files(test_data["path"], tmp_path)
    copy_pseudos("quantum_espresso", tmp_path)
    copy_pseudos("qmcpack", tmp_path)

    for script in test_data["scripts"]:
        script_path = example_root / test_data["path"] / script
        success, message = run_example_script(script_path, test_path)
        assert(success), message

    for code, _filetype, filepath in test_data["files"]:
        success, message = check_generated_files(
            code,
            tmp_path,
            test_data["path"],
            filepath,
            )
        assert(success), message
#end def test_qmcpack_oxygen_dimer


def test_rmg_diamond(tmp_path):
    test_data = dict(
        path = 'rmg/01_diamond_scf',
        scripts = [
            'diamond_scf.py',
            ],
        files = [
            ('rmg', 'input', 'runs/diamond2/scf_gen/scf.in'),
            ('rmg', 'input', 'runs/diamond2/scf_man/scf.in'),
            ],
        )

    test_path = copy_example_files(test_data["path"], tmp_path)
    copy_pseudos("qmcpack", tmp_path)

    for script in test_data["scripts"]:
        script_path = example_root / test_data["path"] / script
        success, message = run_example_script(script_path, test_path)
        assert(success), message

    for code, _filetype, filepath in test_data["files"]:
        success, message = check_generated_files(
            code,
            tmp_path,
            test_data["path"],
            filepath,
            )
        assert(success), message
#end def test_rmg_diamond


def test_pyscf_h2o_ae_hf(tmp_path):
    test_data = dict(
        path = 'pyscf/01_h2o_hf',
        scripts = [
            'h2o_ae_hf.py',
            ],
        files = [
            ('pyscf', 'input', 'runs/h2o_ae_hf/scf.py'),
            ],
        )

    test_path = copy_example_files(
        example_dir   = test_data["path"],
        tmp_dir       = tmp_path,
        template_file = example_root / "pyscf/01_h2o_hf/scf_template.py",
        )

    for script in test_data["scripts"]:
        script_path = example_root / test_data["path"] / script
        success, message = run_example_script(script_path, test_path)
        assert(success), message

    for code, _filetype, filepath in test_data["files"]:
        success, message = check_generated_files(
            code,
            tmp_path,
            test_data["path"],
            filepath,
            )
        assert(success), message
#end def test_pyscf_h2o_ae_hf
