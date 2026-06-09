#! /usr/bin/env python3

from nexus import settings,job,run_project
from nexus import generate_physical_system
from nexus import generate_quantum_package

# note: you must source the QP config file before running this script
#   source /home/ubuntu/apps/qp2/quantum_package.rc

settings(
    results       = '',
    status_only   = 0,
    generate_only = 0,
    sleep         = 3,
    machine       = 'ws16',
    qprc          = '/home/ubuntu/apps/qp2/quantum_package.rc',
    )

scf_job = job(cores=16,threads=16)

system = generate_physical_system(
    structure = 'H2O.xyz',
    )

scf = generate_quantum_package(
    identifier   = 'hf',        # log output goes to hf.out
    path         = 'h2o_ae_hf', # directory to run in
    job          = scf_job,
    system       = system,
    prefix       = 'h2o',       # create/use h2o.ezfio
    run_type     = 'scf',       # qprun scf h2o.ezfio
    ao_basis     = 'cc-pvtz',   # use cc-pvtz basis
    )

run_project()
