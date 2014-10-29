#! /usr/bin/env python

from nexus import settings


settings(
    generate_only = True,                 # only write input files, do not run
    sleep         = 3,                    # check on jobs every 3 seconds
    pseudo_dir    = './pseudopotentials', # path to PP file collection
    machine       = 'node8'               # local machine is an 8 core workstation
    )


settings(
    status_only     = True,                 # only write job status, do not run
    generate_only   = True,                 # only write input files, do not run
    sleep           = 3,                    # check on jobs every 3 seconds
    pseudo_dir      = './pseudopotentials', # path to PP file collection
    local_directory = './'                  # base path for runs and results
    runs            = 'runs',               # runs directory
    results         = 'results',            # results directory
    machine         = 'titan',              # local machine is Titan
    account         = 'MAT123',             # allocation account on Titan
    machine_info    = dict(
        oic5 = dict(
            local_directory = '/home/your_id',
            
            ),
        kraken = dict(
            local_directory = '/nics/b/home/your_id'
            )
        )
    )
