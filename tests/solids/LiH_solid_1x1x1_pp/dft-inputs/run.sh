#!/bin/sh
for kpoint in gamma x arb
do
 pw.x -pw2casino <LiH-${kpoint}-scf.in >& LiH-${kpoint}-scf.out
 pw2qmcpack.x <LiH-${kpoint}-pw2qmcpack.in >& LiH-${kpoint}-pw2qmcpack.out
done
