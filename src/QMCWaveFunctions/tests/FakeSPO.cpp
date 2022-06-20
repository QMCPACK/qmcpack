//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2020 QMCPACK developers.
//
// File developed by: Mark Dewing, mdewing@anl.gov, Argonne National Laboratory
//
// File created by: Mark Dewing, mdewing@anl.gov, Argonne National Laboratory
//////////////////////////////////////////////////////////////////////////////////////


#include "FakeSPO.h"

namespace qmcplusplus
{

FakeSPO::FakeSPO()
{
  className = "FakeSPO";
  a.resize(3, 3);

  a(0, 0) = 2.3;
  a(0, 1) = 4.5;
  a(0, 2) = 2.6;
  a(1, 0) = 0.5;
  a(1, 1) = 8.5;
  a(1, 2) = 3.3;
  a(2, 0) = 1.8;
  a(2, 1) = 4.4;
  a(2, 2) = 4.9;

  v.resize(3);
  v[0] = 1.9;
  v[1] = 2.0;
  v[2] = 3.1;


  a2.resize(4, 4);
  a2(0, 0) = 2.3;
  a2(0, 1) = 4.5;
  a2(0, 2) = 2.6;
  a2(0, 3) = 1.2;
  a2(1, 0) = 0.5;
  a2(1, 1) = 8.5;
  a2(1, 2) = 3.3;
  a2(1, 3) = 0.3;
  a2(2, 0) = 1.8;
  a2(2, 1) = 4.4;
  a2(2, 2) = 4.9;
  a2(2, 3) = 2.8;
  a2(3, 0) = 0.8;
  a2(3, 1) = 4.1;
  a2(3, 2) = 3.2;
  a2(3, 3) = 1.1;

  v2.resize(4, 4);

  v2(0, 0) = 3.2;
  v2(0, 1) = 0.5;
  v2(0, 2) = 5.9;
  v2(0, 3) = 3.7;
  v2(1, 0) = 0.3;
  v2(1, 1) = 1.4;
  v2(1, 2) = 3.9;
  v2(1, 3) = 8.2;
  v2(2, 0) = 3.3;
  v2(2, 1) = 5.4;
  v2(2, 2) = 4.9;
  v2(2, 3) = 2.2;
  v2(3, 1) = 5.4;
  v2(3, 2) = 4.9;
  v2(3, 3) = 2.2;

  gv.resize(4);
  gv[0] = TinyVector<ValueType, DIM>(1.0, 0.0, 0.1);
  gv[1] = TinyVector<ValueType, DIM>(1.0, 2.0, 0.1);
  gv[2] = TinyVector<ValueType, DIM>(2.0, 1.0, 0.1);
  gv[3] = TinyVector<ValueType, DIM>(0.4, 0.3, 0.1);
}

std::unique_ptr<SPOSet> FakeSPO::makeClone() const { return std::make_unique<FakeSPO>(*this); }

void FakeSPO::setOrbitalSetSize(int norbs) { OrbitalSetSize = norbs; }

void FakeSPO::evaluateValue(const ParticleSet& P, int iat, ValueVector& psi)
{
  if (iat < 0)
    for (int i = 0; i < psi.size(); i++)
      psi[i] = 1.2 * i - i * i;
  else if (OrbitalSetSize == 3)
    for (int i = 0; i < 3; i++)
      psi[i] = a(iat, i);
  else if (OrbitalSetSize == 4)
    for (int i = 0; i < 4; i++)
      psi[i] = a2(iat, i);
}

void FakeSPO::evaluateVGL(const ParticleSet& P, int iat, ValueVector& psi, GradVector& dpsi, ValueVector& d2psi)
{
  if (OrbitalSetSize == 3)
  {
    for (int i = 0; i < 3; i++)
    {
      psi[i]  = v[i];
      dpsi[i] = gv[i];
    }
  }
  else if (OrbitalSetSize == 4)
  {
    for (int i = 0; i < 4; i++)
    {
      psi[i]  = v2(iat, i);
      dpsi[i] = gv[i];
    }
  }
}

void FakeSPO::evaluate_notranspose(const ParticleSet& P,
                                   int first,
                                   int last,
                                   ValueMatrix& logdet,
                                   GradMatrix& dlogdet,
                                   ValueMatrix& d2logdet)
{
  if (OrbitalSetSize == 3)
  {
    for (int i = 0; i < 3; i++)
      for (int j = 0; j < 3; j++)
      {
        logdet(j, i)  = a(i, j);
        dlogdet[i][j] = gv[j] + GradType(i);
      }
  }
  else if (OrbitalSetSize == 4)
  {
    for (int i = 0; i < 4; i++)
      for (int j = 0; j < 4; j++)
      {
        logdet(j, i)  = a2(i, j);
        dlogdet[i][j] = gv[j] + GradType(i);
      }
  }
}

} // namespace qmcplusplus
