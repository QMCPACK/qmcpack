#ifndef NLJOB_GPU_H
#define NLJOB_GPU_H

template <typename S>
struct NLjobGPU
{
  int Elec, NumQuadPoints;
  S *R, *QuadPoints, *Ratios;
};

#endif
