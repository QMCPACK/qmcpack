#! /usr/bin/env python3
import os
import sys
import h5py
import numpy as np
import pandas as pd


def read(fdat):
  """ read the scalar.dat file in table format readable by numpy.loadtxt.

   The header line should start with '#' and contain column labels.

  Args:
    dat_fname (str): name of input file
  Return:
    pd.DataFrame: df containing the table of data
  """
  with open(fdat, 'r') as fp:
    header = fp.readline()
  # end with
  cols = header.replace('#', '').split()
  df = pd.read_table(fdat, sep='\s+', comment='#', header=None, names=cols)
  return df
# end def read


def compare_columns_dat_h5(fdat, fh5):
  """ compare mutual data columns in scalar.dat and stat.h5 files

  Args:
    fdat (str): name of scalar.dat file
    fh5  (str): name of stat.h5 file
  Return:
    dict: a dictionary holding mutual columns names as key
  """

  # open database
  df = read(fdat)
  dat_cols = df.columns

  fp = h5py.File(fh5,'r')
  h5_cols = fp.keys()

  # compare mutual columns in .dat v.s. .h5
  agree_map = {}  # keep track of which columns agree
  for col in h5_cols:
    if col not in dat_cols:
      continue

    # check if col agree between .dat and .h5

    # get .h5 values
    h5_loc = os.path.join(col, 'value')
    h5y  = fp[h5_loc][:][:,-1]

    # get .dat values
    daty = df.loc[:,col].values
    agree_map[col] = np.allclose(h5y,daty)
  # end for col
   
  # close database
  fp.close()

  if len(agree_map) == 0:
    raise RuntimeError('%s and %s have no mutual column' % (fdat, fh5))

  return agree_map
# end def


if __name__ == '__main__':

  prefix = sys.argv[1]
  seriesl= [0,1]

  # check Properties
  series_success_map = {}
  for iseries in seriesl:

    # define files to read
    fdat = '%s.s00%d.scalar.dat' % (prefix, iseries)
    fh5  = '%s.s00%d.stat.h5' % (prefix, iseries)

    agree_map = compare_columns_dat_h5(fdat, fh5)
    success = np.all( agree_map.values() )
    series_success_map[iseries] = success
  # end for iseries
  
  all_success = np.all( series_success_map.values() )
  if all_success:
    sys.exit(0)
  else:
    sys.exit(1)
  
# end __main__
