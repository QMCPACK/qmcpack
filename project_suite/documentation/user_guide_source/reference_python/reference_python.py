#! /usr/bin/env python


from project import settings

settings(
    local_directory = './',
    runs            = 'runs',
    results         = 'results',
    pseudo_dir      = './pseudopotentials',
    sleep           = 3,
    generate_only   = False,
    status_only     = False,
    machine         = 'machine_name',
    account         = 'cluster_account',
    machine_info    = dict(
        oic5 = dict(
            local_directory = '/your/path/to/local_directory/on/oic5',
            app_directory   = '/where/you/keep/apps/on/oic5',
            app_directories = {
                'qmcapp'  :'/path/to/qmcapp/on/oic5',
                'pw.x'    :'/path/to/pw.x/on/oic5',
                'myqmcapp':'/path/to/myqmcapp/on/oic5'
                }
            ),
        edison = dict(
            local_directory = '/your/path/to/local_directory/on/edison',
            app_directory   = '/where/you/keep/apps/on/edison',
            app_directories = {
                'qmcapp'  :'/path/to/qmcapp/on/edison',
                'pw.x'    :'/path/to/pw.x/on/edison',
                'myqmcapp':'/path/to/myqmcapp/on/edison'
                }
            )
        )
    )
