import os
import sys
from pathlib import Path

import pytest

from .. import Settings, settings, testing
from ..basisset import BasisSets
from ..developer import obj
from ..gamess import Gamess
from ..generic import NexusUserWarning
from ..machines import Job, Workstation
from ..nexus_base import nexus_config, ShowStatusMode, SimStage
from ..project_manager import ProjectManager
from ..pseudoset import PseudoSet
from ..pwscf import Pwscf
from ..quantum_package import QuantumPackage
from . import NexusTestOrder, isolate_nexus_core

pytestmark = pytest.mark.order(NexusTestOrder.SETTINGS_OPERATION)


@isolate_nexus_core
def test_settings(tmp_path):
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
            'dependent_modes',
            'file_locations',
            'generate_only',
            'graph_sims',
            'indent',
            'load_images',
            'local_directory',
            'monitor',
            'progress_tty',
            'pseudo_dir',
            'quiet',
            'remote_directory',
            'results',
            'runs',
            'skip_submit',
            'sleep',
            'stages',
            'status',
            'timeout',
            'status_only',
            'dynamic',
            'basis_dir',
            'basissets',
            }
        setkeys_check = {
            'command_line',
            'dependent_modes',
            'file_locations',
            'generate_only',
            'graph_sims',
            'indent',
            'load_images',
            'local_directory',
            'monitor',
            'progress_tty',
            'pseudo_dir',
            'quiet',
            'remote_directory',
            'results',
            'runs',
            'skip_submit',
            'sleep',
            'stages',
            'status',
            'timeout',
            'status_only',
            'dynamic',
            'basis_dir',
            'basissets',
            }
        setkeys_allowed = setkeys_check | Settings.allowed_vars

        nckeys  = set(nexus_config.__slots__)
        setkeys = set(settings.keys())

        assert(nckeys==nckeys_check)
        assert(setkeys>=setkeys_check)
        assert(setkeys<=setkeys_allowed)

        for s in nexus_config.__slots__:
            assert(settings[s] == getattr(nexus_config, s))
    #end check_settings_core_noncore

    def check_empty_settings():
        settings(command_line = False)
        settings.command_line     = True
        nexus_config.command_line = True
        check_settings_core_noncore()
        # nexus config has basic run stages and PseudoSet registries are empty
        assert(nexus_config.stages is SimStage.all)

        assert(len(PseudoSet.pseudo_files)==0)
        assert(len(PseudoSet.labeled_pseudosets)==0)
        assert(isinstance(nexus_config.basissets,BasisSets))
        assert(len(nexus_config.basissets)==0)
        nexus_config.restore_defaults()
        assert(nexus_config.basissets is None)
        # other settings objects should be at default also
        aux_defaults()
    #end def_check_empty_settings

    nexus_config.restore_defaults()
    assert(nexus_config.timeout==5*60)
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
    assert(nexus_config.status_only==0)
    assert(nexus_config.generate_only==1)
    assert(nexus_config.timeout==10)
    pseudo_path = str((tmp_path / 'pseudopotentials').resolve())
    assert(nexus_config.pseudo_dir==pseudo_path)
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


@isolate_nexus_core
def test_legacy_and_path_settings(tmp_path):
    local_dir = tmp_path / "local"
    basis_dir = tmp_path / "basis"
    basis_dir.mkdir()

    with pytest.warns(NexusUserWarning, match="verbose"):
        settings(command_line=False, verbose=False)

    settings(
        command_line=False,
        status="standard",
        local_directory=local_dir,
        basis_dir=basis_dir,
        )

    assert(nexus_config.status is ShowStatusMode.all)
    assert(nexus_config.local_directory == str(local_dir))
    assert(nexus_config.basis_dir == str(basis_dir))
    assert(settings.local_directory == str(local_dir))
    assert(settings.basis_dir == str(basis_dir))


@isolate_nexus_core
def test_legacy_settings_preserve_runtime_behavior():
    with pytest.warns(NexusUserWarning, match="verbose"):
        settings(command_line=False, verbose=False)
    assert(nexus_config.quiet)

    with (
        pytest.warns(NexusUserWarning, match="debug"),
        # pytest.warns(NexusUserWarning, match="verbose")
        # Adding the second filter breaks with Pytest 7.4.4
        # See PR #6207
        ):
        settings(command_line=False, verbose=False, debug=True)
    assert(not nexus_config.quiet)

    expected_modes = {
        "none": SimStage(0),
        "setup": SimStage.write_input,
        "send_files": SimStage.send_files,
        "submit": SimStage.submit,
        "get_output": SimStage.get_output,
        "analyze": SimStage.analyze,
        "all": SimStage.all,
        }
    for mode, expected in expected_modes.items():
        with pytest.warns(NexusUserWarning, match="mode"):
            settings(command_line=False, mode=mode)
        assert(nexus_config.stages is expected)

    with pytest.warns(NexusUserWarning, match="mode"):
        settings(command_line=False, mode="submit", stages=["analyze"])
    assert(nexus_config.stages is SimStage.submit)

    settings(command_line=False, stages=[])
    assert(nexus_config.stages is SimStage.all)
    settings(command_line=False, stages="setup")
    assert(nexus_config.stages is SimStage.write_input)

    settings(command_line=False, stages=["submit"], generate_only=True)
    assert(nexus_config.stages is SimStage.submit)
    assert(nexus_config.generate_only)

    settings(command_line=False, sleep=0, timeout=0)
    assert(nexus_config.sleep == 0)
    assert(nexus_config.timeout == 0)
