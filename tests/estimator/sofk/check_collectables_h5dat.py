#! /usr/bin/env python3
import sys
import os
import h5py
import numpy as np
from check_properties_h5dat import read


def get_last_sk(fdat,fh5):
  """ extract S(k) at the longest k vector from scalar.dat and stat.h5

  Args:
    fdat (str): name of scalar.dat file
    fh5  (str): name of stat.h5 file
  Return:
    tuple: (myy, h5y), S(k_max) at each block from scalar.dat and stat.h5
  """

  # get S(k) from scalar.dat
  df = read(fdat)
  sk_cols = [col for col in df.columns if col.startswith('sk')]
  myy = df[sk_cols[-1]].values

  # get S(k) from stat.h5
  fp = h5py.File(fh5, 'r')
  h5y = fp['h5sk/value'][:].T[-1]
  fp.close()

  return myy, h5y
# end def


def show_scalar_trace(data, seriesl):
  import matplotlib.pyplot as plt
  method_map = {0:'VMC',1:'DMC'}
  fig,ax_arr = plt.subplots(1, 2, sharey=True)
  ax_arr[0].set_ylabel('S(k->inf)')
  iplot = 0
  for iseries in seriesl:
    ax = ax_arr[iplot]
    ax.set_title(method_map[iseries])
    ax.set_xlabel('block')
    ax.set_ylim(0.3, 1.2)

    entry = data[iseries]
    daty  = entry['daty']
    h5y   = entry['h5y']

    sline = ax.plot(daty)
    hline = ax.plot(h5y, ls='--', lw=2, alpha=0.8)

    ax.legend(
      handles = [sline[0], hline[0]]
     ,labels  = ['scalar.dat', 'stat.h5']
     ,loc=0
    )

    iplot += 1
  # end for iseries
  plt.show()


if __name__ == '__main__':

  prefix = sys.argv[1]
  seriesl= [0,1]  # a list of series IDs to check

  # check Properties v.s. Collectables
  collectable_success_map = {}
  data = {}
  for iseries in seriesl:

    # define files to read
    fdat = '%s.s00%d.scalar.dat' % (prefix, iseries)
    fh5  = '%s.s00%d.stat.h5' % (prefix, iseries)

    daty, h5y = get_last_sk(fdat, fh5)
    success   = np.allclose(daty, h5y, atol=0.1)
    collectable_success_map[iseries] = success

    # save data for plotting
    data[iseries] = {'daty':daty, 'h5y':h5y}
  # end for 
  all_success = np.all( collectable_success_map.values() )

  if all_success:
    sys.exit(0)
  else:
    #show_scalar_trace(data, seriesl)
    sys.exit(1)

# end __main__
