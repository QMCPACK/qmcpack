
import testing
from testing import divert_nexus,restore_nexus
from testing import value_eq,object_eq


def test_import():
    from nexus import settings,Settings
#end def test_import


def test_settings():
    # test full imports
    import os
    from nexus import settings,Settings,obj
    from nexus_base import nexus_core,nexus_core_defaults
    from nexus_base import nexus_noncore,nexus_noncore_defaults
    from nexus_base import nexus_core_noncore,nexus_core_noncore_defaults
    from pseudopotential import Pseudopotentials
    from basisset import BasisSets
    from machines import Job,Workstation
    from project_manager import ProjectManager
    from gamess import Gamess
    from pwscf import Pwscf
    from quantum_package import QuantumPackage

    testing.check_final_state()

    tpath = testing.setup_unit_test_output_directory('settings','test_settings')

    # divert logging function
    divert_nexus()

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
        nckeys_check = set([
                'command_line','debug', 'dependent_modes', 'emulate',
                'file_locations', 'generate_only', 'graph_sims', 'indent',
                'load_images', 'local_directory', 'mode', 'modes', 'monitor',
                'primary_modes', 'progress_tty', 'pseudo_dir',
                'pseudopotentials', 'remote_directory', 'results', 'runs',
                'skip_submit', 'sleep', 'stages', 'stages_set', 'status',
                'status_modes', 'status_only', 'trace', 'verbose'
                ])
        nnckeys_check = set([
                'basis_dir', 'basissets', 'pseudo_dir', 'pseudopotentials'
                ])
        setkeys_check = set([
                'command_line','basis_dir', 'basissets', 'debug',
                'dependent_modes', 'emulate', 'file_locations', 'generate_only',
                'graph_sims', 'indent', 'load_images', 'local_directory', 'mode',
                'modes', 'monitor', 'primary_modes', 'progress_tty',
                'pseudo_dir', 'pseudopotentials', 'remote_directory', 'results',
                'runs', 'skip_submit', 'sleep', 'stages', 'stages_set', 'status',
                'status_modes', 'status_only', 'trace', 'verbose'
                ])
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
                if isinstance(v1,obj):
                    assert(object_eq(v1,v2))
                else:
                    assert(value_eq(v1,v2))
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
        # nexus core sets basic run stages and Pseudopotentials object
        assert(nexus_core.stages_set==set(nexus_core_defaults.primary_modes))
        assert(isinstance(nexus_core.pseudopotentials,Pseudopotentials))
        assert(len(nexus_core.pseudopotentials)==0)
        nexus_core.stages_set       = set()
        nexus_core.stages           = []
        nexus_core.pseudopotentials = None
        assert(object_eq(nexus_core,nexus_core_defaults))
        # nexus noncore sets Pseudopotentials and BasisSets objects
        assert(isinstance(nexus_noncore.pseudopotentials,Pseudopotentials))
        assert(len(nexus_noncore.pseudopotentials)==0)
        assert(isinstance(nexus_noncore.basissets,BasisSets))
        assert(len(nexus_noncore.basissets)==0)
        nnc_defaults = obj(nexus_noncore_defaults,nexus_core_noncore_defaults)
        nexus_noncore.pseudopotentials = None
        nexus_noncore.basissets        = None
        assert(object_eq(nexus_noncore,nnc_defaults))
        # other settings objects should be at default also
        aux_defaults()
    #end def_check_empty_settings
    
    
    # check that core settings are at default values
    assert(object_eq(nexus_core,nexus_core_defaults))
    assert(object_eq(nexus_noncore,nexus_noncore_defaults))
    assert(object_eq(nexus_core_noncore,nexus_core_noncore_defaults))
    aux_defaults()

    # core settings remain almost at default with empty settings
    check_empty_settings()

    # check that a few basic user settings are applied appropriately
    cwd = os.getcwd()
    os.chdir(tpath)
    dft_pseudos = ['Ni.opt.upf','O.opt.upf']
    qmc_pseudos = ['Ni.opt.xml','O.opt.xml']
    pseudos = dft_pseudos+qmc_pseudos
    pseudo_path = './pseudopotentials'
    if not os.path.exists(pseudo_path):
        os.makedirs(pseudo_path)
        for file in pseudos:
            filepath = os.path.join(pseudo_path,file)
            if not os.path.exists(filepath):
                open(filepath,'w').close()
            #end if
        #end for
    #end if
    settings(
        pseudo_dir    = pseudo_path,
        status_only   = 0,
        generate_only = 1,
        machine       = 'ws16',
        command_line  = False,
        )
    check_settings_core_noncore()
    assert(nexus_core.status_only==0)
    assert(nexus_core.generate_only==1)
    assert(nexus_core.pseudo_dir=='./pseudopotentials')
    assert(len(nexus_core.pseudopotentials)==4)
    assert(set(nexus_core.pseudopotentials.keys())==set(pseudos))
    assert(settings.machine=='ws16')
    assert(Job.machine=='ws16')
    assert(isinstance(ProjectManager.machine,Workstation))
    assert(ProjectManager.machine.name=='ws16')
    os.chdir(cwd)

    # check that a new empty settings works following basic
    check_empty_settings()

    # restore logging function
    restore_nexus()
#end def test_settings
