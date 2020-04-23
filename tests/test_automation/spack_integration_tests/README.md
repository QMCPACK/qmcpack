# Spack Integration Tests

# Description 
An integration test that runs nightly on the Argonne
General Computing Environment (GCE) Linux environment. The script
tests many combinations of variants. It is currently running about 3 times 
a week on `naromero-desktop.cels.anl.gov`

# Usage
Script is designed to run inside a user-level cronjob. `chronic` from more utils is
helpful to silence the script when everything runs correctly. More information
about `chronic` can be found here at https://packages.debian.org/unstable/utils/moreutils

Here is an example of how to add it to your crontab:
```
MAILTO=naromero@anl.gov
0  18 * * 0,3,5 chronic /home/naromero/nightly_anl_gce.sh
```

# Overview of files
## nightly_anl_gce.sh
Deletes all QMCPACK and Quantum Espresso (QE) installations and builds over 200 variants from scratch. 
We also explictly test that the QE patch is applied.

## packages_anl_gce.yaml
To minimize the chance of failure, we include Python and MKL rather than install these with Spack. 
Additionally, for the sake of conciseness in the main test script, we specify concretization and provider preferences 
per the documentation here:
https://spack.readthedocs.io/en/latest/build_settings.html#concretization-preferences

This file should be placed in ` ~/.spack/packages.yaml `

## compilers_anl_gce.yaml
This is the compiler configuration that is very site-specific and is included as an example.

This file should be placed in `~/.spack/linux/compilers.yaml`

If you are using compilers on top of modules, Spack will need to have the correct module name.


## Other Notes
Tested with Spack develop branch with SHA value 25a65638ff19119dcb04db784170eaa44d545a8e
on April 1, 2020.
