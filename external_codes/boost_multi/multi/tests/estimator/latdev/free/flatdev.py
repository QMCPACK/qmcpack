#! /usr/bin/env python3
import h5py
import numpy as np

def corr(trace):
  """ calculate the autocorrelation of a trace of scalar data

  correlation time is defined as the integral of the auto-correlation
   function from t=0 to when the function first reaches 0.

  Args:
    trace (list): should be a 1D iterable array of floating point numbers
  Return:
    float: correlation_time, the autocorrelation time of this trace of scalars
  """
  mu     = np.mean(trace)
  stddev = np.std(trace, ddof=1)
  if np.isclose(stddev, 0):  # easy case
    return np.inf  # infinite correlation for constant trace
  correlation_time = 0.
  for k in range(1, len(trace)):
    # calculate auto_correlation
    auto_correlation = 0.0
    num = len(trace)-k
    auto_correlation = np.dot(trace[:num]-mu, trace[k:]-mu)
    auto_correlation *= 1.0/(num*stddev**2)
    if auto_correlation > 0:
      correlation_time += auto_correlation
    else:
      break
  correlation_time = 1.0 + 2.0*correlation_time
  return correlation_time

def get_latdev_me(fh5, nequil):
  # get raw data
  fp = h5py.File(fh5, 'r')
  data = fp['latdev/value'][()]
  fp.close()

  nb, ncol = data.shape
  npt = nb-nequil
  # calculate mean
  ym = data[nequil:, :].mean(axis=0)
  # calculate standard error
  yc = np.array([corr(data[nequil:, icol]) for icol in range(ncol)])
  neff = npt/yc
  ye = data[nequil:, :].std(ddof=1, axis=0)/neff**0.5
  return ym, ye

def check_latdev(w1, w2, ym, ye, ndim=3):
  expect = np.array([w1**2/2.]*ndim + [w2**2/2.]*ndim)
  zscore = abs(ym-expect)/ye
  fail = np.any(zscore > 3.)
  success = ~fail
  return success, zscore

if __name__ == '__main__':
  # calculate <x^2> <y^2> <z^2>
  nequil = 6
  fh5 = 'mt.s000.stat.h5'
  ym, ye = get_latdev_me(fh5, nequil)

  # check answers
  w1 = 0.3
  w2 = 0.5
  success, zscore = check_latdev(w1, w2, ym, ye)
  if not success:
    print(zscore)

  if success:
    exit(0)
  else:
    exit(1)
# end __main__
