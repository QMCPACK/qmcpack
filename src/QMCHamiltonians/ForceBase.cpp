//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: John R. Gergely,  University of Illinois at Urbana-Champaign
//                    Ken Esler, kpesler@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Raymond Clay III, j.k.rofling@gmail.com, Lawrence Livermore National Laboratory
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: John R. Gergely,  University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////


#include "ForceBase.h"
#include "Particle/DistanceTable.h"
#include "Message/Communicate.h"
#include "Utilities/ProgressReportEngine.h"
#include "Numerics/MatrixOperators.h"
#include "Numerics/DeterminantOperators.h"
#include "OhmmsData/AttributeSet.h"


namespace qmcplusplus
{
ForceBase::ForceBase(ParticleSet& ions, ParticleSet& elns) : FirstForceIndex(-1), tries(0), addionion(true), Ions(ions)
{
  ReportEngine PRE("ForceBase", "ForceBase");
  FirstTime = true;
  Nnuc      = ions.getTotalNum();
  Nel       = elns.getTotalNum();
  //Determines if ion-ion force will be added to electron-ion force in derived force estimators.
  //If false, forces_IonIon=0.0 .
  addionion = true;
  pairName  = elns.getName() + "-" + ions.getName();
  forces.resize(Nnuc);
  forces = 0.0;
  forces_IonIon.resize(Nnuc);
  forces_IonIon = 0.0;
}

void ForceBase::addObservablesF(QMCTraits::PropertySetType& plist)
{
  if (FirstForceIndex < 0)
    FirstForceIndex = plist.size();
  for (int iat = 0; iat < Nnuc; iat++)
  {
    for (int x = 0; x < OHMMS_DIM; x++)
    {
      std::ostringstream obsName;
      obsName << prefix << "_" << iat << "_" << x;
      plist.add(obsName.str());
    }
  }
}

void ForceBase::addObservablesStress(QMCTraits::PropertySetType& plist)
{
  if (FirstForceIndex < 0)
    FirstForceIndex = plist.size();
  for (int i = 0; i < OHMMS_DIM; i++)
    for (int j = i; j < OHMMS_DIM; j++)
    {
      std::ostringstream obsName;
      obsName << prefix << "_" << i << "_" << j;
      plist.add(obsName.str());
    }
}

void ForceBase::registerObservablesF(std::vector<ObservableHelper>& h5list, hid_t gid) const
{
  std::vector<int> ndim(2);
  ndim[0] = Nnuc;
  ndim[1] = OHMMS_DIM;

  h5list.emplace_back(prefix);
  auto& h5o = h5list.back();
  h5o.set_dimensions(ndim, FirstForceIndex);
  h5o.open(gid);
}

void ForceBase::setObservablesF(QMCTraits::PropertySetType& plist)
{
  int index = FirstForceIndex;
  for (int iat = 0; iat < Nnuc; iat++)
  {
    for (int x = 0; x < OHMMS_DIM; x++)
    {
      plist[index] = forces[iat][x];
      index++;
    }
  }
}

void ForceBase::setObservablesStress(QMCTraits::PropertySetType& plist)
{
  int index = FirstForceIndex;
  for (int iat = 0; iat < OHMMS_DIM; iat++)
  {
    for (int jat = iat; jat < OHMMS_DIM; jat++)
    {
      plist[index] = stress(iat, jat);
      index++;
    }
  }
}


void ForceBase::setParticleSetF(QMCTraits::PropertySetType& plist, int offset)
{
  int index = FirstForceIndex + offset;
  for (int iat = 0; iat < Nnuc; iat++)
  {
    for (int x = 0; x < OHMMS_DIM; x++)
    {
      plist[index] = forces[iat][x];
      index++;
    }
  }
}

void ForceBase::setParticleSetStress(QMCTraits::PropertySetType& plist, int offset)
{
  int index = FirstForceIndex + offset;
  for (int iat = 0; iat < OHMMS_DIM; iat++)
  {
    for (int jat = iat; jat < OHMMS_DIM; jat++)
    {
      plist[index] = stress(iat, jat);
      index++;
    }
  }
}

BareForce::BareForce(ParticleSet& ions, ParticleSet& elns) : ForceBase(ions, elns), d_ei_ID(elns.addTable(ions))
{
  name_  = "HF_Force_Base";
  prefix = "HFBase";
}

void BareForce::resetTargetParticleSet(ParticleSet& P) {}

std::unique_ptr<OperatorBase> BareForce::makeClone(ParticleSet& qp, TrialWaveFunction& psi)
{
  return std::make_unique<BareForce>(*this);
}

void BareForce::addObservables(PropertySetType& plist, BufferType& collectables)
{
  addObservablesF(plist);
  my_index_ = FirstForceIndex;
}

BareForce::Return_t BareForce::evaluate(ParticleSet& P)
{
  forces                                    = forces_IonIon;
  const auto& d_ab                          = P.getDistTableAB(d_ei_ID);
  const ParticleSet::Scalar_t* restrict Zat = Ions.Z.first_address();
  const ParticleSet::Scalar_t* restrict Qat = P.Z.first_address();
  //Loop over distinct eln-ion pairs
  for (int jat = 0; jat < d_ab.targets(); jat++)
  {
    const auto& ab_dist  = d_ab.getDistRow(jat);
    const auto& ab_displ = d_ab.getDisplRow(jat);
    for (int iat = 0; iat < d_ab.sources(); iat++)
    {
      real_type rinv = 1.0 / ab_dist[iat];
      real_type r3zz = Qat[jat] * Zat[iat] * rinv * rinv * rinv;
      forces[iat] += r3zz * ab_displ[iat];
    }
  }
  tries++;
  return 0.0;
}

bool BareForce::put(xmlNodePtr cur)
{
  std::string ionionforce("yes");
  OhmmsAttributeSet attr;
  attr.add(prefix, "name");
  attr.add(ionionforce, "addionion");
  attr.put(cur);
  addionion = (ionionforce == "yes" || ionionforce == "true");
  return true;
}

void ForceBase::InitVarReduction(real_type rcut, int _m, int numFuncs)
{
  m    = _m;
  Rcut = rcut;
  std::vector<real_type> h(numFuncs);
  Matrix<real_type> S(numFuncs, numFuncs);
  ck.resize(numFuncs, 0.0);
  real_type R2jp1 = Rcut * Rcut;
  real_type R2m   = 1.0;
  for (int i = 0; i < m; i++)
    R2m *= Rcut;
  for (int j = 1; j <= numFuncs; j++)
  {
    h[j - 1]      = R2jp1 / real_type(j + 1);
    real_type R2k = Rcut;
    for (int k = 1; k <= numFuncs; k++)
    {
      S(k - 1, j - 1) = R2m * R2k * R2jp1 / (real_type)(m + k + j + 1);
      S(k - 1, j - 1) = std::pow(Rcut, (m + k + j + 1)) / (m + k + j + 1.0);
      R2k *= Rcut;
    }
    R2jp1 *= Rcut;
  }
  // fprintf (stderr, "Sij = \n");
  // for (int i=0; i<numFuncs; i++) {
  //   for (int j=0; j<numFuncs; j++)
  // 	fprintf (stderr, " %12.6f ", S(i,j));
  //   fprintf (stderr, "\n");
  // }
  invert_matrix(S, false);
  for (int i = 0; i < numFuncs; i++)
  {
    for (int j = 0; j < numFuncs; j++)
      ck[i] += S(i, j) * h[j];
  }
  FILE* fout = fopen("g_r.dat", "w");
  for (double r = 0.0; r < Rcut; r += 0.001)
    fprintf(fout, "%1.10f %1.10e\n", r, g(r));
  fclose(fout);
  app_log() << "Initialized variance reduction coefs.\n";
}

} // namespace qmcplusplus
