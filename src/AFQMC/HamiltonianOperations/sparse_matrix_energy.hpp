//////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source
// License.  See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by:
// Miguel A. Morales, moralessilva2@llnl.gov
//    Lawrence Livermore National Laboratory
// Alfredo Correa, correaa@llnl.gov
//    Lawrence Livermore National Laboratory
//
// File created by:
// Miguel A. Morales, moralessilva2@llnl.gov
//    Lawrence Livermore National Laboratory
////////////////////////////////////////////////////////////////////////////////

#ifndef  AFQMC_HAMILTONIANOPERATIONS_ENERGY_HPP
#define  AFQMC_HAMILTONIANOPERATIONS_ENERGY_HPP

#include <type_traits>
#include "Configuration.h"
#include "AFQMC/Numerics/ma_operations.hpp"

namespace qmcplusplus
{

namespace afqmc
{

// base == Serial implementation
namespace base
{

/*
 * Calculates the local energy from (already evaluated) mixed density matrices.
 *
 * Vakbl(ak,bl) = sum_i sum_j conj(TrialWfn(i,a)) * conj(TrialWfn(j,b)) * (<ij|kl> - <ij|lk>)
 *    --> (alpha, beta) terms, no exchange term is included (e.g. <ij|lk>)
 * The 2-body contribution to the energy is obtained from:
 * Let G(i,j) {Gmod(i,j)} be the {modified} "Green's Function" of the walker
 *            (defined by the Slater Matrix "W"), then:
 *   G    = conj(TrialWfn) * [transpose(W)*conj(TrialWfn)]^(-1) * transpose(W)
 *   Gmod = [transpose(W)*conj(TrialWfn)]^(-1) * transpose(W)
 *   E_2body = sum_i,j,k,l G(i,k) * (<ij|kl> - <ij|lk>) * G(j,l) + (alpha/beta) + (beta/beta)
 *           = sum a,k,b,l Gmod(a,k) * Vabkl(ak,jl) * Gmod(jl) = Gmod * Vakbl * Gmod
 *   The expression can be evaluated with a sparse matrix-dense vector product,
 *   followed by a dot product between vectors, if we interpret the matrix Gmod(a,k)
 *   as a vector with "linearized" index ak=a*NMO+k.
 */
template< class EMat,
          class MatA,
          class MatB,
          class SpMat
        >
inline void calculate_energy(EMat&& locV, const MatA& Gc, MatB&& Gcloc, const SpMat& Vakbl)
{
  assert(locV.dimensionality==2);
  assert(Gc.shape()[1] == Gcloc.shape()[1]);
  assert(Gc.shape()[0] == Gcloc.shape()[0]);
  assert(Gc.shape()[0] == Vakbl.shape()[0]);
  assert(Gc.shape()[0] == Vakbl.shape()[1]);

  using Type = typename std::decay<MatA>::type::element;
  Type zero = Type(0.);
  Type one = Type(1.);
  Type half = Type(0.5);

  int nwalk = Gc.shape()[1];
  // Vakbl * Gc(bl,nw) = Gcloc(ak,nw)
  ma::product(Vakbl, Gc, std::forward<MatB>(Gcloc));

  for(int i=0, iend=Gc.shape()[0]; i<iend; i++)
    for(int n=0; n<nwalk; n++)
      locV[n][1] += Gc[i][n]*Gcloc[i][n];
  for(int n=0; n<nwalk; n++) locV[n][1] *= half;
}

}

// Shared Memory implementation
namespace shm
{

/*
 * Calculates the local energy from (already evaluated) mixed density matrices.
 *
 * Vakbl(ak,bl) = sum_i sum_j conj(TrialWfn(i,a)) * conj(TrialWfn(j,b)) * (<ij|kl> - <ij|lk>)
 *    --> (alpha, beta) terms, no exchange term is included (e.g. <ij|lk>)
 * The 2-body contribution to the energy is obtained from:
 * Let G(i,j) {Gmod(i,j)} be the {modified} "Green's Function" of the walker
 *            (defined by the Slater Matrix "W"), then:
 *   G    = conj(TrialWfn) * [transpose(W)*conj(TrialWfn)]^(-1) * transpose(W)
 *   Gmod = [transpose(W)*conj(TrialWfn)]^(-1) * transpose(W)
 *   E_2body = sum_i,j,k,l G(i,k) * (<ij|kl> - <ij|lk>) * G(j,l) + (alpha/beta) + (beta/beta)
 *           = sum a,k,b,l Gmod(a,k) * Vabkl(ak,jl) * Gmod(jl) = Gmod * Vakbl * Gmod
 *   The expression can be evaluated with a sparse matrix-dense vector product,
 *   followed by a dot product between vectors, if we interpret the matrix Gmod(a,k)
 *   as a vector with "linearized" index ak=a*NMO+k.
 */
template< class EMat,
          class MatA,
          class MatB,
          class SpMat
        >
inline void calculate_energy(EMat&& locV, const MatA& Gc, MatB&& Gcloc, const SpMat& Vakbl)
{
  // W[nwalk][2][NMO][NAEA]

  assert(locV.dimensionality==2);
  assert(Gc.shape()[1] == Gcloc.shape()[1]);
  assert(Vakbl.shape()[0]  == Gcloc.shape()[0]);
  assert(Gc.shape()[0] == Vakbl.shape()[1]);

  using Type = typename std::decay<MatB>::type::element;
  const Type zero = Type(0.);
  const Type one = Type(1.);
  const Type half = Type(0.5);

  int nwalk = Gc.shape()[1];
  // Vakbl * Gc(bl,nw) = Gcloc(ak,nw)
  ma::product(Vakbl, Gc, std::forward<MatB>(Gcloc));

  // E2(nw) = 0.5*Gc(:,nw)*Gcloc(:,nw)
  int r0 = Vakbl.local_origin()[0];
  for(int i=0, iend=Gcloc.shape()[0]; i<iend; i++)
    for(int n=0; n<nwalk; n++)
      locV[n][1] += Gc[i+r0][n]*Gcloc[i][n];
  for(int n=0; n<nwalk; n++) locV[n][1] *= half;
}

}

}

}

#endif
