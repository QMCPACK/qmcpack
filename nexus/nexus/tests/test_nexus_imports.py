

def test_imports():
    import nexus
    from nexus import settings,job,run_project
    from nexus import ppset
    from nexus import read_structure,read_input
    from nexus import generate_structure
    from nexus import generate_physical_system
    from nexus import generate_qmcpack,generate_convert4qmc
    from nexus import generate_pwscf,generate_pw2qmcpack
    from nexus import generate_gamess
    from nexus import generate_pyscf
    from nexus import generate_quantum_package
    from nexus import generate_vasp
    from nexus import ci,error

    import nexus.versions
    import nexus.nexus_base
    import nexus.testing
    import nexus.execute
    import nexus.memory
    import nexus.generic
    import nexus.developer
    import nexus.unit_converter
    import nexus.periodic_table
    import nexus.numerics
    import nexus.grid_functions
    import nexus.fileio
    import nexus.hdfreader
    import nexus.structure
    import nexus.physical_system
    import nexus.basisset
    import nexus.pseudopotential
    import nexus.machines
    import nexus.simulation
    import nexus.bundle
    import nexus.project_manager
    import nexus.vasp_input
    import nexus.pwscf_input
    import nexus.pwscf_postprocessors
    import nexus.gamess_input
    import nexus.pyscf_input
    import nexus.quantum_package_input
    import nexus.rmg_input
    import nexus.qmcpack_converters
    import nexus.qmcpack_input
    import nexus.vasp_analyzer
    import nexus.pwscf_analyzer
    import nexus.pwscf_postprocessors
    import nexus.gamess_analyzer
    import nexus.pyscf_analyzer
    import nexus.quantum_package_analyzer
    import nexus.rmg_analyzer
    import nexus.qmcpack_converters
    import nexus.qmcpack_analyzer
    import nexus.vasp
    import nexus.pwscf
    import nexus.gamess
    import nexus.pyscf_sim
    import nexus.quantum_package
    import nexus.rmg
    import nexus.pwscf_postprocessors
    import nexus.qmcpack_converters
    import nexus.qmcpack
    import nexus.observables
    
    
    from nexus import settings
    from nexus.machines import job
    from nexus import run_project
    from nexus.pseudopotential import ppset
    from nexus.structure import read_structure
    from nexus import read_input
    from nexus.structure import generate_structure
    from nexus.physical_system import generate_physical_system
    from nexus.qmcpack import generate_qmcpack
    from nexus.qmcpack_converters import generate_convert4qmc
    from nexus.qmcpack_converters import generate_pw2qmcpack
    from nexus.pwscf import generate_pwscf
    from nexus.gamess import generate_gamess
    from nexus.pyscf_sim import generate_pyscf
    from nexus.quantum_package import generate_quantum_package
    from nexus.vasp import generate_vasp
    from nexus.developer import ci,error
    from nexus import read_structure

#end def test_imports
