##################################################################
##  (c) Copyright 2015-  by Jaron T. Krogel                     ##
##################################################################


#====================================================================#
#  nexus_base.py                                                     #
#    Provides base class functionality and access to 'global' data   #
#    for core Nexus classes.                                         #
#                                                                    #
#  Content summary:                                                  #
#    NexusCore                                                       #
#      Base class for core Nexus classes.                            #
#      Data intended to be 'global' among core classes is assigned   #
#      by the Settings class.                                        #
#                                                                    #
#    nexus_core                                                      #
#      Namespace to be accessed by core Nexus classes.               #
#      These are classes that inherit from NexusCore directly.       #
#                                                                    #
#    nexus_noncore                                                   #
#      Namespace to be accessed in read-only fashion by non-core     #
#      classes.                                                      #
#                                                                    #
#====================================================================#
from __future__ import annotations

import os
import pickle
from collections.abc import Collection
from copy import deepcopy
from enum import Flag, auto
from os import PathLike
from pathlib import Path
from pickle import UnpicklingError
from typing import TypeAlias

from .basisset import BasisSets
from .developer import DevBase, log, obj
from .memory import resident
from .nexus_version import nexus_version
from .utilities import path_string

StrPath: TypeAlias = str

class ShowStatusMode(Flag):
    """Flags for controlling the output of a status report.

    See :meth:`~.project_manager.ProjectManager.write_simulation_status`.
    """

    NONE     = auto()
    READY    = auto()
    ACTIVE   = auto()
    FAILED   = auto()
    ALL      = READY | ACTIVE | FAILED
#end class ShowStatusMode


class SimStage(Flag):
    """Flags for controlling which stages of a simulation should be run.

    See :meth:`~.simulation.Simulation.progress`.
    """

    WRITE_INPUT = auto()
    SEND_FILES  = auto()
    SETUP       = WRITE_INPUT | SEND_FILES
    SUBMIT      = auto()
    GET_OUTPUT  = auto()
    ANALYZE     = auto()
    ALL         = SETUP | SUBMIT | GET_OUTPUT | ANALYZE

    @classmethod
    def from_list(cls, items: Collection[str]) -> SimStage:
        if not all(isinstance(item, str) for item in items):
            msg = f"Expected a collection of strs, but got {items}"
            raise TypeError(msg)

        stages = cls(0)
        for stage in items:
            st_up = stage.upper()
            if st_up not in cls.__members__:
                msg = f"Encountered invalid stage '{stage}'"
                raise ValueError(msg)
            stages |= cls[st_up]
        return stages
#end class SimStage


class NexusConfig:
    """Singleton class that represents all of Nexus's configuration settings."""

    __slots__ = (  # noqa: RUF023
        "status_only",
        "status",
        "sleep",
        "timeout",
        "runs",
        "results",
        "local_directory",
        "remote_directory",
        "file_locations",
        "monitor",
        "skip_submit",
        "load_images",
        "stages",
        "dependent_modes",
        "quiet",
        "indent",
        "progress_tty",
        "graph_sims",
        "command_line",
        "dynamic",
        "basis_dir",
        "basissets",
        "pseudo_dir",
        # Legacy
        "modes",
        "mode",
        "generate_only",
        "debug",
        "trace",
        "emulate",
        "primary_modes",
        "status_modes",
    )

    status_only: bool
    """Show status only and exit.

    See :meth:`~.project_manager.ProjectManager.write_simulation_status`.
    """

    status: ShowStatusMode
    """Which jobs to show the status of.

    See :class:`~.ShowStatusMode`, :attr:`~.status_only`.
    """

    sleep: int
    """How long (in seconds) to sleep between polling the queue and checking memory consumption of subprocesses.

    See :meth:`~.project_manager.ProjectManager.run_project`.
    """

    timeout: int
    """Number of seconds to wait for output and error files after a job exits the queue before marking the simulation as failed.

    See :meth:`~.simulation.Simulation.check_status`.
    """

    runs: StrPath
    """Name of the directory that Nexus runs should be placed in.

    See :meth:`~.simulation.Simulation.set_directories`.
    """

    results: StrPath
    """Name of the directory that results from Nexus runs should be copied to.

    If not set (empty string), results will not be stored outside of the runs directory.

    See :meth:`~.simulation.Simulation.set_directories`.
    """

    local_directory: StrPath
    """Directory that Nexus will create files relative to.

    See :meth:`~.simulation.Simulation.set_directories`.
    """

    remote_directory: StrPath
    """Unknown purpose."""

    file_locations: list[StrPath]
    """List of paths to directories that Nexus will look in for files.

    See :meth:`~.simulation.Simulation.send_files`.
    """

    monitor: bool
    """Toggle whether or not Nexus should continue to monitor jobs after submission.

    See :meth:`~.project_manager.ProjectManager.run_project`.
    """

    skip_submit: bool
    """Toggle whether or not to skip actually submitting a job.

    See :meth:`~.simulation.Simulation.submit`.

    Also related to :class:`~.bundle.SimulationBundle`.
    """

    load_images: bool
    """Whether or not to load the simulation image saves to reconstruct the project state.

    See :meth:`~.project_manager.ProjectManager.load_cascades`.
    """

    stages: SimStage
    """Control which simulation stages Nexus should run.

    See :class:`~.SimStage` and :meth:`~.simulation.Simulation.progress`.
    """

    dependent_modes: SimStage
    """Set which stages are required for runtime execution."""

    quiet: bool
    """Disable all Nexus output after initialization."""

    indent: str
    """Indentation base level for Nexus output.

    See :meth:`~.nexus_base.NexusCore.log`.
    """

    progress_tty: bool
    """Toggle printing abbreviated polling messages."""

    graph_sims: bool
    """Optionally create a graph of the simulations that maps their dependency trees.

    See :func:`~nexus.run_project` and :func:`~.simulation.graph_sims`.
    """

    command_line: bool
    """Toggle processing of command line arguments.

    See :meth:`nexus.Settings.__call__`
    """

    dynamic: bool
    """Toggle use of dynamic workflows.

    Used in various places. Some spots to look are in
    :meth:`.simulation.Simulation.__init__`, :func:`~.pwscf.generate_pwscf`,
    :func:`~.qmcpack_converters.generate_pw2qmcpack`, and
    :func:`~.project_manager.workflow_manager`.
    """

    basis_dir: StrPath | None
    """Directory that basis sets are stored in.

    See :attr:`~.basissets`.
    """

    basissets: BasisSets
    """Basis sets found in ``basis_dir`` if it exists.

    See :func:`~.gamess_input.generate_any_gamess_input`.
    """

    pseudo_dir: StrPath | None
    """Directory that pseudopotentials are stored in.

    See :attr:`~.pseudoset.PseudoSet.pseudo_files`,
    :meth:`~.vasp_input.VaspInput.set_potcar`.
    """

    generate_only: bool
    """(LEGACY) Toggle only generating inputs and sending files."""

    def __init__(self):
        self.restore_defaults()

    def restore_defaults(self) -> None:
        self.status_only      = False
        self.status           = ShowStatusMode.NONE
        self.sleep            = 3
        self.timeout          = 5*60
        self.runs             = 'runs'
        self.results          = ''
        self.local_directory  = './'
        self.remote_directory = './'
        self.file_locations   = ['./']
        self.monitor          = True
        self.skip_submit      = False
        self.load_images      = True
        self.stages           = SimStage.ALL
        self.dependent_modes  = SimStage.SUBMIT
        self.quiet            = False
        self.indent           = '  '
        self.progress_tty     = False
        self.graph_sims       = False
        self.command_line     = True
        self.dynamic          = False
        self.basis_dir        = None
        self.basissets        = None
        self.pseudo_dir       = None
        # Legacy
        self.generate_only = False

NEXUS_CONFIG = NexusConfig()


nexus_modules = [mod.stem for mod in Path(__file__).parent.iterdir() if mod.suffix == ".py"]

class NexusUnpickler(pickle.Unpickler):
    """This class is designed for backwards compatibility with pickles generated
    before Nexus was packaged (PR #5700, December 20, 2025).
    It shouldn't touch anything but old Nexus pickles.
    """
    def find_class(self, module, name):
        if module in nexus_modules and "nexus." not in module:
            module = "nexus." + module
        if module == "nexus.generic":
            if name == "obj":
                module = "nexus.developer_tools"
            elif name == "DevBase":
                module = "nexus.developer"

        return super().find_class(module, name)


def write_splash():
    if not hasattr(write_splash, "wrote_splash"):
        splash_text = '''
_____________________________________________________

                     Nexus {}.{}.{}

        (c) Copyright 2012-  Nexus developers

                     Please cite:
  J. T. Krogel Comput. Phys. Commun. 198 154 (2016)
     https://doi.org/10.1016/j.cpc.2015.08.012
_____________________________________________________

'''.format(*nexus_version)
        log(splash_text)
        write_splash.wrote_splash = True
    #end if
#end def write_splash


class NexusCore(DevBase):

    # mutable/dynamic nexus core data
    wrote_something   = False # for pretty printing
    working_directory = None

    def mem_usage(self):
        return int(resident()/1e6)
    #end def mem_usage

    def log(self,*texts,**kwargs):
        """Write output to log file.

        Parameters
        ----------
        *texts
            Strings that will be joined by newlines
        n : int, kwargs
            Spaces to indent
        progress : bool, kwargs
            If ``True`` and output is to a terminal, overwrite and update the
            last line, rather than scrolling.
        """
        if not NEXUS_CONFIG.quiet:
            if len(kwargs)>0:
                n = kwargs['n']
            else:
                n=0
            #end if
            is_progress = kwargs.get('progress',False)
            text=''
            for t in texts:
                text+=str(t)+' '
            #end for
            pad = n*NEXUS_CONFIG.indent
            output_text = pad+text.replace('\n','\n'+pad)
            if NEXUS_CONFIG.progress_tty and is_progress and self._logfile.isatty():
                # spaces to ensure previous line is overwritten.  Need better solution.
                self._logfile.write(output_text+'        \r')
                self._logfile.flush()
            else:
                self._logfile.write(output_text+'\n')
        #end if
        NexusCore.wrote_something = True
    #end def log

    def enter(self, directory: PathLike, *, changedir: bool = True, msg: str = ''):
        """Have Nexus enter a directory and change its current working directory.

        Parameters
        ----------
        directory : PathLike
            Directory to enter. Can be a ``str`` or ``pathlib.Path``
            object.
        changedir : bool, default=True
            Default of ``True`` will change the CWD, setting to ``False``
            will not change the CWD.
        msg : str, optional
            Optional message to pass to the output log.
        """
        NexusCore.working_directory = os.getcwd()
        directory = path_string(directory)

        self.log('    Entering ' + directory, msg)
        if changedir:
            os.chdir(directory)
        #end if
        pad = '      '
        return pad
    #end def enter

    def leave(self):
        os.chdir(NexusCore.working_directory)
    #end def leave

    def load(self, fpath: PathLike | None = None):
        if fpath is None:
            fpath = f'./{type(self).__name__}.p'

        with open(fpath, 'rb') as fobj:
            try:
                tmp = pickle.load(fobj)
            except (ImportError, ModuleNotFoundError):
                fobj.seek(0)
                try:
                    # Old pickles from before Nexus was packaged (PR #5700, December 20 2025)
                    # won't have the correct module path. The custom unpickler will handle this by
                    # prepending "nexus." to the module path
                    tmp = NexusUnpickler(fobj).load()
                except UnpicklingError:
                    # NumPy pickles can use latin1 encoding
                    # They will likely still fail from an underflow since they are not pickle-compliant
                    tmp = NexusUnpickler(fobj).load(encoding='latin1')

        d = self.__dict__
        d.clear()
        for k, v in tmp.__dict__.items():
            d[k] = v
    #end def load
#end class NexusCore


# support dynamic workflows
dynamic_storage = obj(
    simulations         = obj(), # all sims, in dyn proc or not
    simulation_ids      = set(),
    dynamic_processes   = obj(),
    dynamic_process_ids = set(),
    )
