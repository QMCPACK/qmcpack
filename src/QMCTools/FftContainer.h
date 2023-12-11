/////////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois / NCSA Open Source License.
// See LICENSE file in top directory for details .
//
// Copyright ( c ) 2018 QMCPACK developers
//
// File developed by : Luke Shulenburger, lshulen@sandia.gov, Sandia National Laboratories
//
// File created by : Luke Shulenburger, lshulen@sandia.gov, Sandia National Laboratories
/////////////////////////////////////////////////////////////////////////////////////////

#ifndef FFT_CONTAINER_H
#define FFT_CONTAINER_H
#include "fftw3.h"

class FftContainer
{
private:
  fftw_plan plan_;
  int nx_, ny_, nz_;

public:
  int fullSize;
  fftw_complex* rspace;
  fftw_complex* kspace;

public:
  FftContainer(int nx, int ny, int nz);
  double getRsValue(int x, int y, int z, int cplex) const;
  void setRsValue(int x, int y, int z, int cplex, double value);
  double getKsValue(int x, int y, int z, int cplex) const;
  void setKsValue(int x, int y, int z, int cplex, double value);

  int getNx() const { return nx_; }
  int getNy() const { return ny_; }
  int getNz() const { return nz_; }

  int getIndex(int x, int y, int z) const { return z + y * nz_ + x * ny_ * nz_; }
  int getQboxIndex(int x, int y, int z) const { return x + y * nx_ + z * ny_ * nx_; }

  ~FftContainer();
  void executeFFT();
  void fixRsNorm(double factor)
  {
    for (int i = 0; i < fullSize; i++)
    {
      rspace[i][0] *= factor;
      rspace[i][1] *= factor;
    }
  }
  void fixKsNorm(double factor)
  {
    for (int i = 0; i < fullSize; i++)
    {
      kspace[i][0] *= factor;
      kspace[i][1] *= factor;
    }
  }

  double getL2NormRS() const;
  double getL2NormKS() const;
};


FftContainer::FftContainer(int nx, int ny, int nz)
{
  nx_      = nx;
  ny_      = ny;
  nz_      = nz;
  fullSize = nx_ * ny_ * nz_;
  rspace   = (fftw_complex*)fftw_malloc(fullSize * sizeof(fftw_complex));
  kspace   = (fftw_complex*)fftw_malloc(fullSize * sizeof(fftw_complex));
  plan_    = fftw_plan_dft_3d(nx_, ny_, nz_, rspace, kspace, -1, FFTW_ESTIMATE);
}

FftContainer::~FftContainer()
{
  fftw_destroy_plan(plan_);
  fftw_free(rspace);
  fftw_free(kspace);
}

void FftContainer::executeFFT() { fftw_execute(plan_); }

double FftContainer::getL2NormRS() const
{
  double l2norm = 0.0;
  for (int i = 0; i < fullSize; i++)
  {
    l2norm += (rspace[i][0] * rspace[i][0] + rspace[i][1] * rspace[i][1]);
  }
  return l2norm;
}

double FftContainer::getL2NormKS() const
{
  double l2norm = 0.0;
  for (int i = 0; i < fullSize; i++)
  {
    l2norm += (kspace[i][0] * kspace[i][0] + kspace[i][1] * kspace[i][1]);
  }
  return l2norm;
}

double FftContainer::getRsValue(int x, int y, int z, int cplex) const { return rspace[getIndex(x, y, z)][cplex]; }

void FftContainer::setRsValue(int x, int y, int z, int cplex, double value)
{
  rspace[getIndex(x, y, z)][cplex] = value;
}

double FftContainer::getKsValue(int x, int y, int z, int cplex) const { return kspace[getIndex(x, y, z)][cplex]; }

void FftContainer::setKsValue(int x, int y, int z, int cplex, double value)
{
  kspace[getIndex(x, y, z)][cplex] = value;
}


#endif
