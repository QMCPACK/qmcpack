# Reference for the one-rank legacy k-space short test

This reference was generated with the particle-by-particle k-space Jastrow fix
at commit `794459c071f91804ebb0914d2606cf60aebad6f4`, using the complex,
double-precision, CPU-only build on SCITAS Kuma (AMD EPYC 9334, GCC 13.2.0,
OpenMP, `QMC_MPI=OFF`). The input and scalar output are included; only trailing
whitespace in the output header has been removed.

The short test retains 1,000 blocks, 8 steps per block and 16 walkers
(128,000 samples). The reference uses 10,000 blocks (1,280,000 samples), the
same block size, timestep, substeps, warmup and wavefunction, and a fixed random
seed offset of 307 selected before running.

Run from `tests/solids/diamondC_1x1x1_pp` so the existing orbital and
pseudopotential files can be found:

```sh
OMP_NUM_THREADS=16 /path/to/qmcpack_complex qmc-ref/kspace-pbyp/qmc_short_kspace.in.xml
```

Using `simstats` from `tests/scripts/check_scalars.py`, excluding the first two
blocks as in CTest, gives an energy of **-10.4998262734 +/- 0.0005452921 Ha**.
The estimated integrated autocorrelation factor is 1.03671. The expected short
run error scales with the square root of the retained block-count ratio:

```text
error_short = error_reference * sqrt(9998 / 998) = 0.0017259195 Ha
error_check = sqrt(error_reference^2 + error_short^2) = 0.0018100115 Ha
```

The checker uses three times `error_check` as its absolute acceptance bound;
the printed error of each test run does not set that bound. The CMake values
round the reference mean and `error_check` to nine decimal places. Only the
energy reference for `short-diamondC_1x1x1_pp-vmc_sdj_kspace-r1-t16` changes;
its existing variance reference and sample-count check are retained.

An independent real, double-precision run with offset 601 and the same sample
count gives -10.4997128116 +/- 0.0005444510 Ha, agreeing within 0.15 combined
standard errors. Its per-thread seed range does not overlap the reference's.
Both runs were planned in advance; neither was retried to obtain agreement.
