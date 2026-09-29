Nexus QE+QMCPACK Example 5: DMC restart
========================================

These scripts restart a single-twist diamond DMC calculation through
``dependencies=(qmc1,'restart')``.  The initial calculation must enable
checkpointing.

``diamond_lda_dmc_restart_same_dir.py`` places both QMCPACK simulations in
one directory and continues the original QMCPACK project series.
``diamond_lda_dmc_restart_separate_dirs.py`` uses separate directories and
references the initial checkpoint by a relative path.
