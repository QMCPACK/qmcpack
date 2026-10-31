import json
import os
import shutil
import sys
from pathlib import Path
from subprocess import PIPE, Popen

import numpy as np
import pytest

from . import TEST_DIR

pytestmark = pytest.mark.real_run

NEXUS_ROOT = TEST_DIR.parent.parent # qmcpack/nexus
EXAMPLE_ROOT  = NEXUS_ROOT / "nexus/examples"
TEST_ROOT     = NEXUS_ROOT / "nexus/tests"
REFERENCE_DIR = TEST_ROOT / "reference/dummy_workflows"

def run_example_script(script: Path, test_path: Path) -> tuple[bool, str]:
    script = script.resolve() # Absolute path helps in debugging
    old_cwd = Path.cwd()
    os.chdir(test_path)

    script_command = (
        f"PYTHONPATH={NEXUS_ROOT} "
        f"{sys.executable} "
        f"{script} "
        "--sleep=0.01 "
        "--disable-logging "
        f"--pseudo-dir={EXAMPLE_ROOT / 'qmcpack/pseudopotentials'}"
    )
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
    tmp_dir: Path,
    filepath: str,
    ) -> None:
    __tracebackhide__ = True

    gen_file = tmp_dir / filepath
    if not gen_file.exists():
        msg = f"Input file was not generated: {gen_file}"
        raise FileNotFoundError(msg)
#end def check_generated_files


@pytest.mark.skipif(
    shutil.which("pw.x") is None,
    reason="Must have pw.x in your $PATH!"
)
def test_qe_no_dependencies(tmp_path):
    script_path = REFERENCE_DIR / "qe_no_dependencies.py"
    success, message = run_example_script(script_path, tmp_path)
    assert(success), message

    required_files = (
        "analysis.json",
        "pwscf_output",
        "scf.err",
        "scf.in",
        "scf.out",
        "sim_scf",
    )

    runs_dir = tmp_path / "runs/scf"
    for file in required_files:
        check_generated_files(runs_dir, file)

    ref = {
        "energy": -9.084907703102045,
        "kpoints": [[0.0, 0.0, 0.0]],
        "kweights": [2.0],
        "eigenvalues": [
            [
                -7.285414577886064,
                19.00581872844667,
                19.006207160066,
                19.007033680918152
            ]
        ],
        "occupations": [[1.0, 1.0, 1.0, 1.0]],
        "Ef": 19.007033680918152,
        "Evbm": None,
        "Ecbm": None,
        "band_gap": None,
        "fractional_occs": False
    }
    calc_analysis = runs_dir / "analysis.json"
    with open(calc_analysis, "r") as out:
        calc = json.load(out)

    assert(calc["Evbm"] is ref["Evbm"])
    assert(calc["Ecbm"] is ref["Ecbm"])
    assert(calc["band_gap"] is ref["band_gap"])
    assert(calc["fractional_occs"] is ref["fractional_occs"])

    np.testing.assert_allclose(calc["occupations"], ref["occupations"])
    np.testing.assert_allclose(calc["kpoints"], ref["kpoints"])
    np.testing.assert_allclose(calc["kweights"], ref["kweights"])
    np.testing.assert_allclose(calc["eigenvalues"], ref["eigenvalues"])

    np.testing.assert_allclose(calc["energy"], ref["energy"])
    np.testing.assert_allclose(calc["Ef"], ref["Ef"])
#end def test_qe_no_dependencies

@pytest.mark.skipif(
    shutil.which("pw.x") is None,
    reason="Must have pw.x in your $PATH!"
)
def test_qe_with_dependencies(tmp_path):
    script_path = REFERENCE_DIR / "qe_with_dependencies.py"
    success, message = run_example_script(script_path, tmp_path)
    assert(success), message

    scf_required_files = (
        "scf_analysis.json",
        "scf_output",
        "scf.err",
        "scf.in",
        "scf.out",
        "sim_scf",
    )
    relax_required_files = (
        "relax_analysis.json",
        "relax_output",
        "relax.err",
        "relax.in",
        "relax.out",
        "sim_relax",
    )

    scf_dir = tmp_path / "runs/scf"
    for file in scf_required_files:
        check_generated_files(scf_dir, file)

    relax_dir = tmp_path / "runs/relax"
    for file in relax_required_files:
        check_generated_files(relax_dir, file)

    scf_ref = {
        "energy": -13.82229685442133,
        "kpoints": [[0.0, 0.0, 0.0]],
        "kweights": [2.0],
        "eigenvalues": [
            [
                -20.138878830563087,
                14.206752904785214,
                16.709686167962005,
                19.247823996746305
            ]
        ],
        "occupations": [[1.0, 1.0, 1.0, 1.0]],
        "Ef": 19.247823996746305,
        "Evbm": None,
        "Ecbm": None,
        "band_gap": None,
        "fractional_occs": False,
        "initial_structure": {
            "axes": [
                [2.000000002237739, 0.0, 0.0],
                [0.0, 2.000000002237739, 0.0],
                [0.0, 0.0, 2.000000002237739]
            ],
            "pos": [
                [ 0.868386705058762,   0.8466856980787831, 1.0000000037647554],
                [ 1.7840978760513462,  0.1772073710430425, 1.0000000037647554],
                [-0.04009316412654221, 1.4953092254088007, 1.0000000037647554]
            ],
            "elem": ["O", "H", "H"]
        }
    }

    scf_analysis_file = scf_dir / "scf_analysis.json"
    with open(scf_analysis_file, "r") as out:
        scf_calc = json.load(out)

    assert(scf_calc["Evbm"] is scf_ref["Evbm"])
    assert(scf_calc["Ecbm"] is scf_ref["Ecbm"])
    assert(scf_calc["band_gap"] is scf_ref["band_gap"])
    assert(scf_calc["fractional_occs"] is scf_ref["fractional_occs"])

    np.testing.assert_allclose(scf_calc["occupations"], scf_ref["occupations"])
    np.testing.assert_allclose(scf_calc["kpoints"], scf_ref["kpoints"])
    np.testing.assert_allclose(scf_calc["kweights"], scf_ref["kweights"])
    np.testing.assert_allclose(scf_calc["eigenvalues"], scf_ref["eigenvalues"])

    np.testing.assert_allclose(scf_calc["energy"], scf_ref["energy"])
    np.testing.assert_allclose(scf_calc["Ef"], scf_ref["Ef"])

    scf_calc_structure = scf_calc["initial_structure"]
    scf_ref_structure = scf_ref["initial_structure"]
    np.testing.assert_allclose(scf_calc_structure["axes"], scf_ref_structure["axes"])
    np.testing.assert_allclose(scf_calc_structure["pos"], scf_ref_structure["pos"])
    assert(scf_calc_structure["elem"] ==  scf_ref_structure["elem"])

    relax_ref = {
        "energy": -13.82522899985655,
        "kpoints": [[0.0, 0.0, 0.0]],
        "kweights": [2.0],
        "eigenvalues": [
            [
                -20.50277734756552,
                14.095760420126403,
                16.816093567788148,
                19.235170317369718
            ]
        ],
        "occupations": [[1.0, 1.0, 1.0, 1.0]],
        "Ef": 19.235170317369718,
        "Evbm": None,
        "Ecbm": None,
        "band_gap": None,
        "fractional_occs": False,
        "initial_structure": {
            "axes": [
                [2.000000002237739, 0.0, 0.0],
                [0.0, 2.000000002237739, 0.0],
                [0.0, 0.0, 2.000000002237739]
            ],
            "pos": [
                [0.6123914147458271, 0.5192022975846599, 1.0000000037647554],
                [1.5963914175189946, 0.5192022975846599, 1.0000000037647554],
                [0.40360858471874417, 1.480797704653079, 1.0000000037647554]
            ],
            "elem": ["O", "H", "H"]
        },
        "relaxed_structure": {
            "axes": [
                [2.000000002237739, 0.0, 0.0],
                [0.0, 2.000000002237739, 0.0],
                [0.0, 0.0, 2.000000002237739]
            ],
            "pos": [
                [ 0.8683867074979552,   0.8466856989756067,  1.0000000037647554],
                [ 1.7840978744713785,   0.17720737284288535, 1.0000000037647554],
                [-0.040093164985765324, 1.4953092280039055,  1.0000000037647554]
            ],
            "elem": ["O", "H", "H"]
        }
    }

    relax_analysis_file = relax_dir / "relax_analysis.json"
    with open(relax_analysis_file, "r") as out:
        relax_calc = json.load(out)

    assert(relax_calc["Evbm"] is relax_ref["Evbm"])
    assert(relax_calc["Ecbm"] is relax_ref["Ecbm"])
    assert(relax_calc["band_gap"] is relax_ref["band_gap"])
    assert(relax_calc["fractional_occs"] is relax_ref["fractional_occs"])

    np.testing.assert_allclose(relax_calc["occupations"], relax_ref["occupations"])
    np.testing.assert_allclose(relax_calc["kpoints"], relax_ref["kpoints"])
    np.testing.assert_allclose(relax_calc["kweights"], relax_ref["kweights"])
    np.testing.assert_allclose(relax_calc["eigenvalues"], relax_ref["eigenvalues"])

    np.testing.assert_allclose(relax_calc["energy"], relax_ref["energy"])
    np.testing.assert_allclose(relax_calc["Ef"], relax_ref["Ef"])

    relax_calc_structure_initial = relax_calc["initial_structure"]
    relax_ref_structure_initial = relax_ref["initial_structure"]
    np.testing.assert_allclose(
        relax_calc_structure_initial["axes"], relax_ref_structure_initial["axes"]
    )
    np.testing.assert_allclose(
        relax_calc_structure_initial["pos"], relax_ref_structure_initial["pos"]
    )
    assert(relax_calc_structure_initial["elem"] ==  relax_ref_structure_initial["elem"])

    relax_calc_structure_relaxed = relax_calc["relaxed_structure"]
    relax_ref_structure_relaxed = relax_ref["relaxed_structure"]
    np.testing.assert_allclose(
        relax_calc_structure_relaxed["axes"], relax_ref_structure_relaxed["axes"]
    )
    np.testing.assert_allclose(
        relax_calc_structure_relaxed["pos"], relax_ref_structure_relaxed["pos"]
    )
    assert(relax_calc_structure_relaxed["elem"] ==  relax_ref_structure_relaxed["elem"])


    # Check that the relaxed structure was used for the SCF
    np.testing.assert_allclose(
        scf_calc_structure["axes"], relax_calc_structure_relaxed["axes"]
    )
    np.testing.assert_allclose(
        scf_calc_structure["pos"], relax_calc_structure_relaxed["pos"]
    )
    assert(scf_calc_structure["elem"] ==  relax_calc_structure_relaxed["elem"])
#end def test_qe_no_dependencies
