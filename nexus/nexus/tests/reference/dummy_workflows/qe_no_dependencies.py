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

structure = generate_structure(lattice="diamond", cell="prim")
system = generate_physical_system(
    structure=structure,
    **pseudos["espresso"].get_Zeff({"C"})
)

scf = generate_pwscf(
    identifier   = 'scf',
    path         = 'scf',
    job          = job(cores=1,app='pw.x'),
    input_type   = 'generic',
    calculation  = 'scf',
    input_dft    = 'lda',
    ecutwfc      = 10,
    conv_thr     = 1e-3,
    system       = system,
    pseudos      = pseudos,
    kgrid        = (1,1,1),
    kshift       = (0,0,0),
)

run_project()

analyzer = analyze_output(scf)
analyzer_data = {
    "energy": analyzer.energy(),
    "kpoints": analyzer.kpoints().tolist(),
    "kweights": analyzer.kweights().tolist(),
    "eigenvalues": analyzer.eigenvalues().tolist(),
    "occupations": analyzer.occupations().tolist(),
    "Ef": analyzer.Ef(),
    "Evbm": analyzer.Evbm(),
    "Ecbm": analyzer.Ecbm(),
    "band_gap": analyzer.band_gap(),
    "fractional_occs": analyzer.fractional_occs(),
}

with open(Path(scf.locdir) / "analysis.json", "w") as out:
    json.dump(analyzer_data, out, indent=4)
