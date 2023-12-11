#! /usr/bin/env python3

from pyscf.pbc import df, scf

$system

gdf = df.FFTDF(cell,kpts)
gdf.auxbasis = 'weigend'

mf = scf.KRHF(cell,kpts).density_fit()
mf.exxdiv  = 'ewald'
mf.with_df = gdf
mf.kernel()
