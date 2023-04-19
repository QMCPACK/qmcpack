#!/bin/bash

#TODO: add tests to do all of this with cmake/ctest

# dependencies:
#   pyscf, numpy, h5py
#   ${QMCROOT}/src/QMCTools in path and pythonpath (PyscfSphToCart.py in path, PyscfToQmcpack.py in pythonpath)
#   convert4qmc, qmcpack, qmca in path

# generate he_${l}_sph.h5
python3 he2_1x1x1_sph2cart.py > pyscf.out

export OMP_NUM_THREADS=1
for l in sp sd #sf sg sh si
do
  # generate he_${l}_cart.h5
  PyscfSphToCart.py -o he_${l}_cart.h5 he_${l}_sph.h5 > s2c.${l}.out
  mkdir $l
  cd $l
  for bas in cart sph
  do
    jobname=he_${l}_${bas}
    cp ../${jobname}.h5 .
    convert4qmc -orbitals ${jobname}.h5 -nojastrow > ${jobname}.convert.out
    sed -e "s/SYSTEMID/${jobname}/g" ../ref.vmc.xml > ${jobname}.vmc.xml
    qmcpack ${jobname}.vmc.xml > ${jobname}.out
    qmca -e0 -q e --fp=16.8e ${jobname}.s000.scalar.dat | awk -F= '{print $2}' > ${jobname}.e.out
  done

  if $(cmp --silent he_${l}_{sph,cart}.e.out)
  then
    echo "he_${l} sph-to-cart: pass"
  else
    echo "he_${l} sph-to-cart: fail"
  fi
  cd ..
done


