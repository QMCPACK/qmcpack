from inspect import signature
from enum import IntEnum, auto
from pathlib import Path
from copy import deepcopy
import functools
from nexus.nexus_base import NEXUS_CONFIG
from nexus.generic import generic_settings
from nexus.pseudoset import PseudoSet
from nexus.simulation import Simulation

# qmcpack/nexus/nexus/tests/
TEST_DIR = Path(__file__).resolve().parent

class FakeLog:
    def __init__(self):
        self.reset()

    def reset(self):
        self.s = ""

    def write(self,s):
        self.s += s

    def close(self):
        pass

    def contents(self):
        return self.s


def divert_nexus_log():
    """Create a fake logging object to divert Nexus's output."""
    logging_storage = {
        'devlog': generic_settings.devlog,
        }
    logfile = FakeLog()
    generic_settings.devlog   = logfile
    return logfile, logging_storage


def restore_nexus_log(logging_storage: dict):
    """Restore Nexus's logging to the state stored in ``logging_storage``."""
    generic_settings.devlog   = logging_storage.pop('devlog')


def isolate_nexus_core(test_func = None):
    """Isolate changes in ``NEXUS_CONFIG`` for a test function."""

    needs_tmp_path = "tmp_path" in str(signature(test_func))

    @functools.wraps(test_func)
    def wrap_path(tmp_path):
        pseudo_files = deepcopy(PseudoSet.pseudo_files)
        labeled_pseudosets = deepcopy(PseudoSet.labeled_pseudosets)
        logfile, logging_storage = divert_nexus_log()
        try:
            test_func(tmp_path)
        finally:
            NEXUS_CONFIG.restore_defaults()
            PseudoSet.pseudo_files = pseudo_files
            PseudoSet.labeled_pseudosets = labeled_pseudosets
            restore_nexus_log(logging_storage)
            Simulation.clear_all_sims()

    @functools.wraps(test_func)
    def wrap():
        pseudo_files = deepcopy(PseudoSet.pseudo_files)
        labeled_pseudosets = deepcopy(PseudoSet.labeled_pseudosets)
        logfile, logging_storage = divert_nexus_log()
        try:
            test_func()
        finally:
            NEXUS_CONFIG.restore_defaults()
            PseudoSet.pseudo_files = pseudo_files
            PseudoSet.labeled_pseudosets = labeled_pseudosets
            restore_nexus_log(logging_storage)
            Simulation.clear_all_sims()

    if needs_tmp_path:
        return wrap_path
    else:
        return wrap


def create_pseudo_files(
    tmp_dir: Path,
    pseudos: list[str],
    pseudo_strs: list[str | None] | None = None
    ):
    """Create pseudopotential files and register them with PseudoSet.

    This function must be called in a function that has been decorated
    with ``@isolate_nexus_core``.

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

    for pseudo, text in zip(pseudos, pseudo_strs):
        pseudo_file = pseudo_dir / pseudo
        pseudo_file.write_text(text)


    PseudoSet.pseudo_files = {
        pseudo.name:str(pseudo.resolve()) for pseudo in pseudo_dir.iterdir()
        if pseudo.is_file()
        }
    PseudoSet.labeled_pseudosets = {}
    NEXUS_CONFIG.pseudo_dir    = str(pseudo_dir)


def register_pseudo_files(pseudos: list[str]):
    """Register synthetic pseudopotential paths for input-generation tests."""
    PseudoSet.pseudo_files.update({
        pseudo:str(Path(pseudo).resolve()) for pseudo in pseudos
        })
#end def register_pseudo_files


class NexusTestOrder(IntEnum):
    """Test order for Nexus testing.
    
    This dictates the order that the tests are run in, reflecting the
    inheritance hierarchy that Nexus has, so the first tests to fail are
    going to be indicative of where the actual root problem is.
    """

    TESTING                         = auto()
    EXECUTE                         = auto()
    MEMORY                          = auto()
    UTILITIES                       = auto()
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
    PSEUDOSET                       = auto()
    PSEUDOPOTENTIAL                 = auto()
    NEXUS_BASE                      = auto()
    ERROR_KEYS                      = auto()
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
