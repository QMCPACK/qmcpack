#ifndef FFT_CONTAINER_H
#define FFT_CONTAINER_H
#include "fftw3.h"

class fftContainer {
 private:
  fftw_plan Plan;
  int nx, ny, nz;
 public:
  int fullSize;
  fftw_complex *rspace;
  fftw_complex *kspace;
 public:
  fftContainer(int _nx, int _ny, int _nz);
  double getRsValue(int x, int y, int z, int cplex) const;
  void setRsValue(int x, int y, int z, int cplex, double value);
  double getKsValue(int x, int y, int z, int cplex) const;
  void setKsValue(int x, int y, int z, int cplex, double value);

  int getNx() const { return nx; }
  int getNy() const { return ny; }
  int getNz() const { return nz; }

  int getIndex(int x, int y, int z) const {
    return z + y * nz + x * ny*nz;
  }
  int getQboxIndex(int x, int y, int z) const {
    return x + y * nx + z * ny*nx;
  }

  ~fftContainer();
  void executeFFT();
  void fixRsNorm(double factor) {
    for (int i = 0; i < fullSize; i++) {
      rspace[i][0] *= factor;
      rspace[i][1] *= factor;
    }
  }
  void fixKsNorm(double factor) {
    for (int i = 0; i < fullSize; i++) {
      kspace[i][0] *= factor;
      kspace[i][1] *= factor;
    }
  }

  double getL2NormRS() const;
  double getL2NormKS() const;
};


fftContainer::fftContainer(int _nx, int _ny, int _nz) {
  nx = _nx;
  ny = _ny;
  nz = _nz;
  fullSize = nx*ny*nz;
  rspace = (fftw_complex*) fftw_malloc(fullSize * sizeof(fftw_complex));
  kspace = (fftw_complex*) fftw_malloc(fullSize * sizeof(fftw_complex));
  Plan = fftw_plan_dft_3d(nx, ny, nz, rspace, kspace, -1, FFTW_ESTIMATE);
}

fftContainer::~fftContainer() {
  fftw_destroy_plan(Plan);
  fftw_free(rspace);
  fftw_free(kspace);
}

void fftContainer::executeFFT() {
  fftw_execute(Plan);
}

double fftContainer::getL2NormRS() const {
  double l2norm = 0.0;
  for (int i = 0; i < fullSize; i++) {
    l2norm += (rspace[i][0]*rspace[i][0] + rspace[i][1]*rspace[i][1]);
  }
  return l2norm;
}

double fftContainer::getL2NormKS() const {
  double l2norm = 0.0;
  for (int i = 0; i < fullSize; i++) {
    l2norm += (kspace[i][0]*kspace[i][0] + kspace[i][1]*kspace[i][1]);
  }
  return l2norm;
}

double fftContainer::getRsValue(int x, int y, int z, int cplex) const {
  return rspace[getIndex(x,y,z)][cplex];
}

void fftContainer::setRsValue(int x, int y, int z, int cplex, double value) {
  rspace[getIndex(x,y,z)][cplex] = value;
}

double fftContainer::getKsValue(int x, int y, int z, int cplex) const {
  return kspace[getIndex(x,y,z)][cplex];
}

void fftContainer::setKsValue(int x, int y, int z, int cplex, double value) {
  kspace[getIndex(x,y,z)][cplex] = value;
}


#endif
