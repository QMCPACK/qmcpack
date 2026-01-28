#! /usr/bin/env python3

def find_nexus_modules():
    import os
    import sys
    # Prepend the assumed path of a folder containing the closely coupled nexus module to module search path
    # It points to the top Nexus directory not necessarily the top QMCPACK directory.
    nexus_lib = os.path.realpath(os.path.join(__file__,'..','..'))
    sys.path.insert(0,nexus_lib)

    # check import of nexus modules
    try:
        import nexus
    except (ImportError, ModuleNotFoundError):
        raise ImportError(
            "Nexus module is required but cannot be imported!"
        )
    nexus_lib2 = os.path.realpath(os.path.join(nexus.__file__,'..','..'))

    # ensure closely coupled nexus module is in use.
    if nexus_lib != nexus_lib2:
        exe_full_path = os.path.abspath(sys.argv[0])
        raise Exception('Broken Nexus installation! Nexus Python scripts are required to reside in the directory that contains the Nexus module.'
                        + '\nNexus Python script         : ' + exe_full_path
                        + '\nNexus module root directory : ' + nexus_lib2)

#end def find_nexus_modules
