#ifndef CUDA_SPLINE_H
#define CUDA_SPLINE_H

#include "QMCWaveFunctions/Jastrow/BsplineFunctor.h"

namespace qmcplusplus
{

template<typename T>
struct CudaSpline
{
  gpu::device_vector<T> coefs;
  T rMax;

  template<typename T2>
  void set (BsplineFunctor<T2> &func)
  {
    int num_coefs = func.SplineCoefs.size();
    gpu::host_vector<T> coefs_h(num_coefs);
    for (int i=0; i<num_coefs; i++)
    {
      coefs_h[i] = func.SplineCoefs[i];
      // app_log() << "coefs_h[" << i << "] = " << coefs_h[i] << endl;
    }
    coefs = coefs_h;
    rMax = func.cutoff_radius;
  }

  void check_coefs()
  {
    gpu::host_vector<T> coefs_h(coefs.size());
    coefs_h = coefs;
    for (int i=0; i<coefs.size(); i++)
      if (isinf(coefs_h[i]))
        fprintf (stderr, "coefs[%d] = %1.5f\n", i, coefs_h[i]);
  }

  template<typename T2>
  CudaSpline (BsplineFunctor<T2> &func) : coefs("CudaSpline::coefs")
  {
    set (func);
  }
};
}

#endif
