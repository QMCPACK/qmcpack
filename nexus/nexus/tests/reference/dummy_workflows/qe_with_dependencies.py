import json
from pathlib import Path

from nexus import (
    analyze_output,
    generate_physical_system,
    generate_pseudoset,
    generate_pwscf,
    generate_structure,
    job,
    run_project,
    settings,
)

settings(
    disable_logging=True,
    machine='ws',
    sleep=0.01,
)

pseudos = generate_pseudoset(
    # pseudo_dir provided on command line
    extension={"rmg": ".xml"},
)

structure = generate_structure(
    type       = "trimer",
    trimer     = ["O", "H", "H"],
    units      = "A",
    separation = [0.984, 0.984],
    angle      = 102.25,
    Lbox       = 2.0,
)
structure.center_molecule()

mol = generate_physical_system(
    structure  = structure,
    **pseudos["espresso"].get_Zeff({"H", "O"})
)

relax = generate_pwscf(
    identifier      = "relax",
    path            = "relax",
    job             = job(cores=1, app="pw.x"),
    input_type      = "generic",
    calculation     = "relax",
    prefix          = "relax",
    outdir          = "relax_output",
    etot_conv_thr   = 1.00e-1,
    forc_conv_thr   = 1.00e-1,
    ecutwfc         = 10,
    input_dft       = "pbe",
    assume_isolated = "mp",
    conv_thr        = 1e-1,
    ion_dynamics    = "bfgs",
    system          = mol,
    pseudos         = pseudos,
    kgrid           = (1,1,1),
    kshift          = (0,0,0),
)

scf = generate_pwscf(
    identifier      = "scf",
    path            = "scf",
    job             = job(cores=1, app="pw.x"),
    input_type      = "generic",
    dependencies    = (relax, "structure"),
    calculation     = "scf",
    prefix          = "scf",
    outdir          = "scf_output",
    ecutwfc         = 10,
    input_dft       = "pbe",
    assume_isolated = "mp",
    conv_thr        = 1e-1,
    system          = mol,
    pseudos         = pseudos,
    kgrid           = (1,1,1),
    kshift          = (0,0,0),
)

run_project()

relax_analyzer = analyze_output(relax)
initial_structure = relax_analyzer.initial_structure()
relaxed_structure = relax_analyzer.relaxed_structure()

relax_analyzer_data = {
    "energy": relax_analyzer.energy(),
    "kpoints": relax_analyzer.kpoints().tolist(),
    "kweights": relax_analyzer.kweights().tolist(),
    "eigenvalues": relax_analyzer.eigenvalues().tolist(),
    "occupations": relax_analyzer.occupations().tolist(),
    "Ef": relax_analyzer.Ef(),
    "Evbm": relax_analyzer.Evbm(),
    "Ecbm": relax_analyzer.Ecbm(),
    "band_gap": relax_analyzer.band_gap(),
    "fractional_occs": relax_analyzer.fractional_occs(),
    "initial_structure": {
        "axes": initial_structure.axes.tolist(),
        "pos": initial_structure.pos.tolist(),
        "elem": initial_structure.elem.tolist(),
    },
    "relaxed_structure": {
        "axes": relaxed_structure.axes.tolist(),
        "pos": relaxed_structure.pos.tolist(),
        "elem": relaxed_structure.elem.tolist(),
    },
}

with open(Path(relax.locdir) / "relax_analysis.json", "w") as out:
    json.dump(relax_analyzer_data, out, indent=4)


scf_analyzer = analyze_output(scf)
scf_initial_structure = scf_analyzer.initial_structure()

scf_analyzer_data = {
    "energy": scf_analyzer.energy(),
    "kpoints": scf_analyzer.kpoints().tolist(),
    "kweights": scf_analyzer.kweights().tolist(),
    "eigenvalues": scf_analyzer.eigenvalues().tolist(),
    "occupations": scf_analyzer.occupations().tolist(),
    "Ef": scf_analyzer.Ef(),
    "Evbm": scf_analyzer.Evbm(),
    "Ecbm": scf_analyzer.Ecbm(),
    "band_gap": scf_analyzer.band_gap(),
    "fractional_occs": scf_analyzer.fractional_occs(),
    "initial_structure": {
        "axes": scf_initial_structure.axes.tolist(),
        "pos": scf_initial_structure.pos.tolist(),
        "elem": scf_initial_structure.elem.tolist(),
    },
}

with open(Path(scf.locdir) / "scf_analysis.json", "w") as out:
    json.dump(scf_analyzer_data, out, indent=4)
