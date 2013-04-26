#ifndef CUDA_WALKER_H
#define CUDA_WALKER_H

#include <vector>
#include "multi_bspline_cuda_s.h"

class cuda_determinant
{
public:
  int N;
  float *A, *Atran, *Ainv;
  float *Ainv_delta, *Ainv_colk;
  float *new_row, *delta;

  void resize(int N);
  cuda_determinant(int N);
  cuda_determinant();
};


class cuda_walker
{
public:
  int N[2];
  float *R;
  cuda_determinant dets[2];
  bool accept;
  void resize(int nup, int ndown);
};


class cuda_population
{
private:
  const int MaxPop;
  float **A_list_d, **Ainv_list_d, **delta_list_d;
  float **Ainv_delta_list_d, **Ainv_colk_list_d;
  float *ratios_d;
  float *pos_d;
  std::vector<float*> A_vec, Ainv_vec, delta_vec,
      Ainv_delta_vec, Ainv_colk_vec;
  std::vector<float> ratio_vec, pos_vec;
  vector<cuda_walker> walkers;
  // Number of up and down electrons
  int num_elecs[2];
  multi_UBspline_3d_s_cuda *multi_spline;
public:
  void calc_new_row (int elec);
  void calc_ratios (int elec);
  void update_determinants(int elec);

  cuda_population();
};

#endif
