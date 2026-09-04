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
from ..nexus_base import NEXUS_CONFIG, SimStage, ShowStatusMode

@isolate_nexus_core
def test_settings(tmp_path):
    # test full imports
    import os
    from nexus import settings,Settings
    from ..developer import DevBase
    
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

    def check_empty_settings():
        settings(command_line = False)
        settings.command_line     = True
        NEXUS_CONFIG.command_line = True
        # nexus config has basic run stages and PseudoSet registries are empty
        assert(NEXUS_CONFIG.stages is SimStage.ALL)

        assert(len(PseudoSet.pseudo_files)==0)
        assert(len(PseudoSet.labeled_pseudosets)==0)
        assert(isinstance(NEXUS_CONFIG.basissets,BasisSets))
        assert(len(NEXUS_CONFIG.basissets)==0)
        NEXUS_CONFIG.restore_defaults()
        assert(NEXUS_CONFIG.basissets is None)
        # other settings objects should be at default also
        aux_defaults()
    #end def_check_empty_settings

    NEXUS_CONFIG.restore_defaults()
    assert(NEXUS_CONFIG.timeout==5*60)
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
    assert(NEXUS_CONFIG.status_only==0)
    assert(NEXUS_CONFIG.generate_only==1)
    assert(NEXUS_CONFIG.timeout==10)
    pseudo_path = str((tmp_path / 'pseudopotentials').resolve())
    assert(NEXUS_CONFIG.pseudo_dir==pseudo_path)
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
