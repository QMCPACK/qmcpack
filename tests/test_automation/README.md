Two things reside in this directory:

# scheduled nightly and weekly test builds
These run tests of QMCPACK on various machines
automatically. Users should run ctest directly, as described in the
manual.

These scripts are designed to be called from cron to execute the
various nightly, weekly tests of QMCPACK. They are custom to
particular machines, file systems, and users. 

A *copy* of these scripts should be used for automation. At every
update check for unsafe operations ("rm *") or script actions that can
go awry e.g. if a filesystem is full or unavailable.

Note that if you are copying/reusing these files, they might not
correspond to exactly those used in production, since these are
sometimes tweaked for runtime, minor software version updates etc.

# CI configuration
Most of these files begin with jenkins. Scrutinize any updates to them
carefully there is basically no reason they should be touched except by maintainers.
They depend on a common jenkins pipeline library pulled from QMCPACK/qmcpack_shared_jenkins.
Updates to that library don't seem to propagate faster than ~10 minutes.
