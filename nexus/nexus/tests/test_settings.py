import pytest
import sys
from . import NexusTestOrder
pytestmark = pytest.mark.order(NexusTestOrder.SETTINGS_OPERATION)

from .. import settings
from ..developer import obj

from pathlib import Path
from . import isolate_nexus_core
from .. import testing
from ..testing import object_eq


@isolate_nexus_core
def test_settings(tmp_path):
    # test full imports
    import os
    from nexus import settings,Settings
    from ..developer import DevBase
    from ..nexus_base import nexus_core,nexus_core_defaults
    from ..nexus_base import nexus_noncore,nexus_noncore_defaults
    from ..nexus_base import nexus_core_noncore,nexus_core_noncore_defaults
    from ..pseudoset import PseudoSet
    from ..basisset import BasisSets
    from ..machines import Job,Workstation
    from ..project_manager import ProjectManager
    from ..gamess import Gamess
    from ..pwscf import Pwscf
    from ..quantum_package import QuantumPackage

    testing.check_final_state()

    def aux_defaults():
        # check that Job and ProjectManager settings are at default values
        assert(Job.machine is None)
        assert(ProjectManager.machine is None)

        # check that Gamess, Pwscf, and Quantum Package settings are at default values
        assert(Gamess.ericfmt is None)
        assert(Gamess.mcppath is None)
        assert(Pwscf.vdw_table is None)
        assert(QuantumPackage.qprc is None)
    #end def aux_defaults

    def check_settings_core_noncore():
        nckeys_check = {
                'command_line',
                'file_locations', 'generate_only', 'graph_sims', 'indent',
                'load_images', 'local_directory', 'monitor',
                'progress_tty', 'pseudo_dir',
                'remote_directory', 'results', 'runs',
                'skip_submit', 'sleep', 'timeout',
                'status_only', 'dynamic'
                }
        nnckeys_check = {
                'basis_dir', 'basissets', 'pseudo_dir'
                }
        setkeys_check = {
                'command_line', 'basis_dir', 'basissets',
                'file_locations', 'generate_only',
                'graph_sims', 'indent', 'load_images', 'local_directory',
                'monitor', 'progress_tty',
                'pseudo_dir', 'remote_directory', 'results',
                'runs', 'skip_submit', 'sleep',
                'timeout', 'status_only', 'dynamic'
                }
        setkeys_allowed = setkeys_check | Settings.allowed_vars

        nckeys  = set(nexus_core.keys())
        nnckeys = set(nexus_noncore.keys())
        setkeys = set(settings.keys())

        assert(nckeys==nckeys_check)
        assert(nnckeys==nnckeys_check)
        assert(setkeys>=setkeys_check)
        assert(setkeys<=setkeys_allowed)

        pairs = [(settings,nexus_core),
                 (settings,nexus_noncore),
                 (nexus_core,nexus_noncore)
                 ]
        for o1,o2 in pairs:
            shared_keys = set(o1.keys()) & set(o2.keys())
            for k in shared_keys:
                v1 = o1[k]
                v2 = o2[k]
                if isinstance(v1,(obj,DevBase)):
                    assert(object_eq(v1,v2))
                else:
                    assert(v1 == v2)
                #end if
            #end for
        #end for
    #end check_settings_core_noncore

    def check_empty_settings():
        settings(
            command_line = False,
            )
        settings.command_line   = True
        nexus_core.command_line = True
        check_settings_core_noncore()
        # PseudoSet registries are empty
        assert(len(PseudoSet.pseudo_files)==0)
        assert(len(PseudoSet.labeled_pseudosets)==0)
        assert(object_eq(nexus_core,nexus_core_defaults))
        # nexus noncore sets a BasisSets object
        assert(isinstance(nexus_noncore.basissets,BasisSets))
        assert(len(nexus_noncore.basissets)==0)
        nnc_defaults = obj(**nexus_noncore_defaults)
        nnc_defaults.update(**nexus_core_noncore_defaults)
        nexus_noncore.basissets        = None
        assert(object_eq(nexus_noncore,nnc_defaults))
        # other settings objects should be at default also
        aux_defaults()
    #end def_check_empty_settings
    
    
    # check that core settings are at default values
    assert(object_eq(nexus_core,nexus_core_defaults))
    assert(nexus_core.timeout==5*60)
    assert(object_eq(nexus_noncore,nexus_noncore_defaults))
    assert(object_eq(nexus_core_noncore,nexus_core_noncore_defaults))
    aux_defaults()

    # core settings remain almost at default with empty settings
    check_empty_settings()

    # check that a few basic user settings are applied appropriately
    cwd = Path.cwd()
    os.chdir(tmp_path)
    dft_pseudos = ['Ni.opt.upf','O.opt.upf']
    qmc_pseudos = ['Ni.opt.xml','O.opt.xml']
    pseudos = dft_pseudos+qmc_pseudos
    pseudo_path = './pseudopotentials'
    if not os.path.exists(pseudo_path):
        os.makedirs(pseudo_path)
        for file in pseudos:
            filepath = Path(pseudo_path) / file
            if not filepath.exists():
                filepath.touch()
            #end if
        #end for
    #end if
    settings(
        pseudo_dir    = pseudo_path,
        status_only   = 0,
        generate_only = 1,
        timeout       = 10,
        machine       = 'ws16',
        command_line  = False,
        )
    check_settings_core_noncore()
    assert(nexus_core.status_only==0)
    assert(nexus_core.generate_only==1)
    assert(nexus_core.timeout==10)
    pseudo_path = str((tmp_path / 'pseudopotentials').resolve())
    assert(nexus_core.pseudo_dir==pseudo_path)
    assert(PseudoSet.pseudo_files=={
        pseudo:str((Path(pseudo_path)/pseudo).resolve()) for pseudo in pseudos
        })
    assert(len(PseudoSet.labeled_pseudosets)==0)
    assert(settings.machine=='ws16')
    assert(Job.machine=='ws16')
    assert(isinstance(ProjectManager.machine,Workstation))
    assert(ProjectManager.machine.name=='ws16')
    os.chdir(cwd)

    # check that a new empty settings works following basic
    check_empty_settings()
#end def test_settings


@isolate_nexus_core
def test_command_line_timeout():
    argv = sys.argv
    try:
        sys.argv = ['nexus_script.py','--timeout=12.5']
        script_settings = obj()
        settings.process_command_line_settings(script_settings)
    finally:
        sys.argv = argv
    #end try

    assert(script_settings.timeout==12.5)
#end def test_command_line_timeout
