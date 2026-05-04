from inspect import signature
from enum import IntEnum, auto
from pathlib import Path
from copy import deepcopy
import functools
from nexus.nexus_base import (
    nexus_core,
    nexus_noncore,
    nexus_core_noncore,
    nexus_noncore_defaults,
)
from nexus.generic import generic_settings, object_interface
from nexus.pseudopotential import Pseudopotentials
from nexus.simulation import Simulation

# qmcpack/nexus/nexus/tests/
TEST_DIR = Path(__file__).resolve().parent

NEXUS_CORE_KEYS = (
    "local_directory",
    "remote_directory",
    "mode",
    "stages",
    "stages_set",
    "status",
    "sleep",
    "file_locations",
    "pseudo_dir",
    "pseudopotentials",
    "runs",
    "results",
)
NEXUS_NONCORE_KEYS = (
    "pseudo_dir",
    "pseudopotentials",
)

def divert_nexus_core():
    """Store Nexus's core and noncore keys and return them."""
    nexus_core_storage = {}
    for key in NEXUS_CORE_KEYS:
        nexus_core_storage[key] = nexus_core[key]
        nexus_core[key] = deepcopy(nexus_core[key])

    nexus_noncore_storage = {}
    for key in NEXUS_NONCORE_KEYS:
        if key in nexus_noncore:
            nexus_noncore_storage[key] = nexus_noncore[key]
            nexus_noncore[key] = deepcopy(nexus_noncore[key])

    return nexus_core_storage, nexus_noncore_storage


def restore_nexus_core(nexus_core_storage: dict, nexus_noncore_storage: dict):
    """Use the keys in ``nexus_core_storage`` and ``nexus_noncore_storage`` to restore state."""

    for key in NEXUS_CORE_KEYS:
        nexus_core[key] = nexus_core_storage.pop(key)

    for key in NEXUS_NONCORE_KEYS:
        if key in nexus_noncore_storage:
            nexus_noncore[key] = nexus_noncore_storage.pop(key)
        elif key in nexus_noncore:
            del nexus_noncore[key]

    for key in list(nexus_noncore.keys()):
        if key not in nexus_noncore_defaults:
            del nexus_noncore[key]

    nexus_core_noncore.pseudopotentials = None

    assert len(nexus_noncore_storage) == 0, "Nexus Core keys have not been properly reset!"
    assert len(nexus_core_storage) == 0,    "Nexus NonCore keys have not been properly reset!"


class FakeLog:
    def __init__(self):
        self.reset()

    def reset(self):
        self.s = ""

    def write(self,s):
        self.s += s

    def close(self):
        None

    def contents(self):
        return self.s


def divert_nexus_log():
    """Create a fake logging object to divert Nexus's output."""
    logging_storage = {
        'devlog': generic_settings.devlog,
        'objlog': object_interface._logfile,
    }
    logfile = FakeLog()
    generic_settings.devlog   = logfile
    object_interface._logfile = logfile
    return logfile, logging_storage


def restore_nexus_log(logging_storage: dict):
    """Restore Nexus's logging to the state stored in ``logging_storage``."""
    generic_settings.devlog   = logging_storage.pop('devlog')
    object_interface._logfile = logging_storage.pop('objlog')


def isolate_nexus_core(test_func = None):
    """Isolate changes in ``nexus_core`` for a test function."""

    needs_tmp_path = "tmp_path" in str(signature(test_func))

    @functools.wraps(test_func)
    def wrap_path(tmp_path):
        nexus_core_storage, nexus_noncore_storage = divert_nexus_core()
        logfile, logging_storage = divert_nexus_log()
        try:
            test_func(tmp_path)
            test_err = None
        except Exception as err:
            test_err = err

        restore_nexus_core(nexus_core_storage, nexus_noncore_storage)
        restore_nexus_log(logging_storage)
        Simulation.clear_all_sims()
        if test_err is not None:
            raise test_err

    @functools.wraps(test_func)
    def wrap():
        nexus_core_storage, nexus_noncore_storage = divert_nexus_core()
        logfile, logging_storage = divert_nexus_log()
        try:
            test_func()
            test_err = None
        except Exception as err:
            test_err = err

        restore_nexus_core(nexus_core_storage, nexus_noncore_storage)
        restore_nexus_log(logging_storage)
        Simulation.clear_all_sims()
        if test_err is not None:
            raise test_err

    if needs_tmp_path:
        return wrap_path
    else:
        return wrap


def create_pseudo_files(
    tmp_dir: Path,
    pseudos: list[str],
    pseudo_strs: list[str | None] | None = None
):
    """Create pseudopotential files and add them to the global pseudopotentials.

    This function must be called in a function that has been decorated
    with ``@isolate_nexus_core(needs_tmp_path=True)``.

    Parameters
    ----------
    tmp_dir : Path
        Path to the temporary directory of a test.
    pseudos : list of str
        List of pseudopotential name(s). These are created at the path
        ``tmp_dir/pseudopotential/<pseudos>``.
    pseudo_strs : list of str or None, optional
        Text to write into the pseudopotential file(s). Must have the
        same length as ``pseudos``.
    """

    if pseudo_strs is None:
        pseudo_strs = [""] * len(pseudos)
    elif len(pseudo_strs) != len(pseudos):
        raise ValueError(
            "Test improperly written!\n"
            "Must have pseudo text or `None` for every pseudo file!"
        )

    pseudo_dir = tmp_dir / "pseudopotentials"
    pseudo_dir.mkdir(parents=True)

    new_pseudos = []
    for pseudo, text in zip(pseudos, pseudo_strs):
        pseudo_file = pseudo_dir / pseudo
        pseudo_file.write_text(text)
        new_pseudos.append(pseudo_file)


    pseudopotentials = Pseudopotentials(new_pseudos)
    nexus_core.pseudopotentials    = pseudopotentials
    nexus_noncore.pseudopotentials = pseudopotentials
    nexus_core.pseudo_dir    = str(pseudo_dir)
    nexus_noncore.pseudo_dir = str(pseudo_dir)


class NexusTestOrder(IntEnum):
    """Test order for Nexus testing.
    
    This dictates the order that the tests are run in, reflecting the
    inheritance hierarchy that Nexus has, so the first tests to fail are
    going to be indicative of where the actual root problem is.
    """

    VERSIONS                        = auto() # 1
    REQUIRED_DEPENDENCIES           = auto() # 2
    OPTIONAL_DEPENDENCIES           = auto() # etc.
    NEXUS_IMPORTS                   = auto()
    TESTING                         = auto()
    EXECUTE                         = auto()
    MEMORY                          = auto()
    GENERIC_OPERATION               = auto()
    DEVELOPER                       = auto()
    UNIT_CONVERTER                  = auto()
    PERIODIC_TABLE                  = auto()
    NUMERICS                        = auto()
    GRID_FUNCTIONS                  = auto()
    FILEIO                          = auto()
    HDFREADER                       = auto()
    XMLREADER                       = auto()
    STRUCTURE                       = auto()
    PHYSICAL_SYSTEM                 = auto()
    BASISSET                        = auto()
    PSEUDOPOTENTIAL                 = auto()
    NEXUS_BASE                      = auto()
    MACHINES                        = auto()
    SIMULATION                      = auto()
    BUNDLE                          = auto()
    PROJECT_MANAGER                 = auto()
    SETTINGS_OPERATION              = auto()
    VASP_INPUT                      = auto()
    PWSCF_INPUT                     = auto()
    PWSCF_POSTPROCESSOR_INPUT       = auto()
    GAMESS_INPUT                    = auto()
    PYSCF_INPUT                     = auto()
    QUANTUM_PACKAGE_INPUT           = auto()
    RMG_INPUT                       = auto()
    QMCPACK_CONVERTER_INPUT         = auto()
    QMCPACK_INPUT                   = auto()
    VASP_ANALYZER                   = auto()
    PWSCF_ANALYZER                  = auto()
    PWSCF_POSTPROCESSOR_ANALYZERS   = auto()
    GAMESS_ANALYZER                 = auto()
    PYSCF_ANALYZER                  = auto()
    QUANTUM_PACKAGE_ANALYZER        = auto()
    RMG_ANALYZER                    = auto()
    QMCPACK_CONVERTER_ANALYZERS     = auto()
    QMCPACK_ANALYZER                = auto()
    VASP_SIMULATION                 = auto()
    PWSCF_SIMULATION                = auto()
    GAMESS_SIMULATION               = auto()
    PYSCF_SIMULATION                = auto()
    QUANTUM_PACKAGE_SIMULATION      = auto()
    RMG_SIMULATION                  = auto()
    PWSCF_POSTPROCESSOR_SIMULATIONS = auto()
    QMCPACK_CONVERTER_SIMULATIONS   = auto()
    QMCPACK_SIMULATION              = auto()
    OBSERVABLES                     = auto()
    NXS_REDO                        = auto()
    NXS_SIM                         = auto()
    QMC_FIT                         = auto()
    QDENS                           = auto()
    QDENS_RADIAL                    = auto()
    QMCA                            = auto()
    USER_EXAMPLES                   = auto()
