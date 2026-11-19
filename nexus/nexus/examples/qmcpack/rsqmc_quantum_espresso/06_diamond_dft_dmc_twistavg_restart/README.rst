Nexus QE+QMCPACK Example 6: Twist-averaged DMC restart
======================================================

These scripts extend Example 5 to a 2x2x2 twist grid.  Nexus matches each
restart checkpoint to the corresponding downstream twist.

``diamond_lda_dmc_twistavg_restart_same_dir.py`` continues each QMCPACK
project series in one directory.  The ``separate_dirs`` variant places the
initial and restarted calculations in different directories.
