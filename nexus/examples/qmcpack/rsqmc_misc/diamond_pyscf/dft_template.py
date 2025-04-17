#! /usr/bin/env python3

from pyscf.pbc import df, scf

$system

gdf = df.GDF(cell,kpts)
gdf.auxbasis = 'weigend'
gdf.build()

mf = scf.KRKS(cell,kpts).density_fit()
mf.xc      ='b3lyp'
mf.tol     = 1e-10
mf.exxdiv  = 'ewald'
mf.with_df = gdf
mf.kernel()
