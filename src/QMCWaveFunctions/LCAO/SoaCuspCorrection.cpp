//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2021 QMCPACK developers.
//
// File developed by: Mark Dewing, mdewing@anl.gov, Argonne National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@intel.com, Intel Corp.
//////////////////////////////////////////////////////////////////////////////////////


/** @file SoaCuspCorrection.cpp
 */
#include "SoaCuspCorrection.h"
#include "SoaCuspCorrectionBasisSet.h"

namespace qmcplusplus
{
SoaCuspCorrection::SoaCuspCorrection(ParticleSet& ions, ParticleSet& els) : myTableIndex(els.addTable(ions))
{
  NumCenters = ions.getTotalNum();
  NumTargets = els.getTotalNum();
  LOBasisSet.resize(NumCenters);
}

SoaCuspCorrection::SoaCuspCorrection(const SoaCuspCorrection& a) = default;

void SoaCuspCorrection::setBasisSetSize(int nbs)
{
  BasisSetSize = nbs;
  //THIS NEEDS TO BE FIXE for OpenMP
  myVGL.resize(5, BasisSetSize);
}

inline void SoaCuspCorrection::evaluateVGL(const ParticleSet& P, int iat, VGLVector_t& vgl)
{
  myVGL = 0.0;

  const auto& d_table = P.getDistTableAB(myTableIndex);
  const auto& dist    = (P.activePtcl == iat) ? d_table.getTempDists() : d_table.getDistRow(iat);
  const auto& displ   = (P.activePtcl == iat) ? d_table.getTempDispls() : d_table.getDisplRow(iat);
  for (int c = 0; c < NumCenters; c++)
  {
    if (LOBasisSet[c])
    {
      LOBasisSet[c]->evaluate_vgl(dist[c], displ[c], myVGL[0], myVGL[1], myVGL[2], myVGL[3], myVGL[4]);
    }
  }

  {
    const auto v_in  = myVGL[0];
    const auto gx_in = myVGL[1];
    const auto gy_in = myVGL[2];
    const auto gz_in = myVGL[3];
    const auto l_in  = myVGL[4];
    auto v_out       = vgl.data(0);
    auto gx_out      = vgl.data(1);
    auto gy_out      = vgl.data(2);
    auto gz_out      = vgl.data(3);
    auto l_out       = vgl.data(4);
    for (size_t i = 0; i < BasisSetSize; ++i)
    {
      v_out[i] += v_in[i];
      gx_out[i] += gx_in[i];
      gy_out[i] += gy_in[i];
      gz_out[i] += gz_in[i];
      l_out[i] += l_in[i];
    }
  }
}

void SoaCuspCorrection::evaluate_vgl(const ParticleSet& P,
                                     int iat,
                                     ValueVector_t& psi,
                                     GradVector_t& dpsi,
                                     ValueVector_t& d2psi)
{
  myVGL = 0.0;

  const auto& d_table = P.getDistTableAB(myTableIndex);
  const auto& dist    = (P.activePtcl == iat) ? d_table.getTempDists() : d_table.getDistRow(iat);
  const auto& displ   = (P.activePtcl == iat) ? d_table.getTempDispls() : d_table.getDisplRow(iat);
  for (int c = 0; c < NumCenters; c++)
  {
    if (LOBasisSet[c])
    {
      LOBasisSet[c]->evaluate_vgl(dist[c], displ[c], myVGL[0], myVGL[1], myVGL[2], myVGL[3], myVGL[4]);
    }
  }

  const auto v_in  = myVGL[0];
  const auto gx_in = myVGL[1];
  const auto gy_in = myVGL[2];
  const auto gz_in = myVGL[3];
  const auto l_in  = myVGL[4];
  for (size_t i = 0; i < BasisSetSize; ++i)
  {
    psi[i] += v_in[i];
    dpsi[i][0] += gx_in[i];
    dpsi[i][1] += gy_in[i];
    dpsi[i][2] += gz_in[i];
    d2psi[i] += l_in[i];
  }
}

void SoaCuspCorrection::evaluate_vgl(const ParticleSet& P,
                                     int iat,
                                     int idx,
                                     ValueMatrix_t& psi,
                                     GradMatrix_t& dpsi,
                                     ValueMatrix_t& d2psi)
{
  myVGL = 0.0;

  const auto& d_table = P.getDistTableAB(myTableIndex);
  const auto& dist    = (P.activePtcl == iat) ? d_table.getTempDists() : d_table.getDistRow(iat);
  const auto& displ   = (P.activePtcl == iat) ? d_table.getTempDispls() : d_table.getDisplRow(iat);
  for (int c = 0; c < NumCenters; c++)
  {
    if (LOBasisSet[c])
    {
      LOBasisSet[c]->evaluate_vgl(dist[c], displ[c], myVGL[0], myVGL[1], myVGL[2], myVGL[3], myVGL[4]);
    }
  }

  const auto v_in  = myVGL[0];
  const auto gx_in = myVGL[1];
  const auto gy_in = myVGL[2];
  const auto gz_in = myVGL[3];
  const auto l_in  = myVGL[4];
  for (size_t i = 0; i < BasisSetSize; ++i)
  {
    psi[idx][i] += v_in[i];
    dpsi[idx][i][0] += gx_in[i];
    dpsi[idx][i][1] += gy_in[i];
    dpsi[idx][i][2] += gz_in[i];
    d2psi[idx][i] += l_in[i];
  }
}

void SoaCuspCorrection::evaluateV(const ParticleSet& P, int iat, ValueType* restrict vals)
{
  ValueType* tmp_vals = myVGL[0];

  std::fill_n(tmp_vals, myVGL.size(), 0.0);

  const auto& d_table = P.getDistTableAB(myTableIndex);
  const auto& dist    = (P.activePtcl == iat) ? d_table.getTempDists() : d_table.getDistRow(iat);

  //THIS IS SERIAL, only way to avoid this is to use myVGL
  for (int c = 0; c < NumCenters; c++)
  {
    if (LOBasisSet[c])
    {
      LOBasisSet[c]->evaluate(dist[c], tmp_vals);
    }
  }

  { //collect
    const auto v_in = myVGL[0];
    for (size_t i = 0; i < BasisSetSize; ++i)
    {
      vals[i] += v_in[i];
    }
  }
}

void SoaCuspCorrection::add(int icenter, std::unique_ptr<COT> aos) { LOBasisSet[icenter].reset(aos.release()); }

} // namespace qmcplusplus
