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
#include "Particle/ParticleSet.h"
#include "Message/Communicate.h"
#include "Utilities/ProgressReportEngine.h"
#include "Numerics/MatrixOperators.h"
#include "Numerics/DeterminantOperators.h"

namespace qmcplusplus
{
using Real = ForceBase::Real;

ForceBase::ForceBase(ParticleSet& ions, ParticleSet& elns)
    : first_force_index_(-1),
      n_nuc_(ions.getTotalNum()),
      n_el_(elns.getTotalNum()),
      tries_(0),
      first_time_(true),
      add_ion_ion_(true),
      ions_(ions)
{
  ReportEngine PRE("ForceBase", "ForceBase");
  pair_name_ = elns.getName() + "-" + ions.getName();
  forces_.resize(n_nuc_);
  forces_ = 0.0;
  forces_ion_ion_.resize(n_nuc_);
  forces_ion_ion_ = 0.0;
}

ForceBase::~ForceBase() {}

void ForceBase::addObservablesF(QMCTraits::PropertySetType& plist)
{
  if (first_force_index_ < 0)
    first_force_index_ = plist.size();
  for (int iat = 0; iat < n_nuc_; iat++)
  {
    for (int x = 0; x < OHMMS_DIM; x++)
    {
      std::ostringstream obsName;
      obsName << prefix_ << "_" << iat << "_" << x;
      plist.add(obsName.str());
    }
  }
}

void ForceBase::addObservablesStress(QMCTraits::PropertySetType& plist)
{
  if (first_force_index_ < 0)
    first_force_index_ = plist.size();
  for (int i = 0; i < OHMMS_DIM; i++)
    for (int j = i; j < OHMMS_DIM; j++)
    {
      std::ostringstream obsName;
      obsName << prefix_ << "_" << i << "_" << j;
      plist.add(obsName.str());
    }
}

void ForceBase::registerObservablesF(std::vector<ObservableHelper>& h5list, hdf_archive& file) const
{
  std::vector<int> ndim(2);
  ndim[0] = n_nuc_;
  ndim[1] = OHMMS_DIM;

  h5list.emplace_back(hdf_path{prefix_});
  auto& h5o = h5list.back();
  h5o.set_dimensions(ndim, first_force_index_);
}

void ForceBase::setObservablesF(QMCTraits::PropertySetType& plist)
{
  int index = first_force_index_;
  for (int iat = 0; iat < n_nuc_; iat++)
  {
    for (int x = 0; x < OHMMS_DIM; x++)
    {
      plist[index] = forces_[iat][x];
      index++;
    }
  }
}

void ForceBase::setObservablesStress(QMCTraits::PropertySetType& plist)
{
  int index = first_force_index_;
  for (int iat = 0; iat < OHMMS_DIM; iat++)
  {
    for (int jat = iat; jat < OHMMS_DIM; jat++)
    {
      plist[index] = stress_(iat, jat);
      index++;
    }
  }
}


void ForceBase::setParticleSetF(QMCTraits::PropertySetType& plist, int offset)
{
  int index = first_force_index_ + offset;
  for (int iat = 0; iat < n_nuc_; iat++)
  {
    for (int x = 0; x < OHMMS_DIM; x++)
    {
      plist[index] = forces_[iat][x];
      index++;
    }
  }
}

void ForceBase::setParticleSetStress(QMCTraits::PropertySetType& plist, int offset)
{
  int index = first_force_index_ + offset;
  for (int iat = 0; iat < OHMMS_DIM; iat++)
  {
    for (int jat = iat; jat < OHMMS_DIM; jat++)
    {
      plist[index] = stress_(iat, jat);
      index++;
    }
  }
}

void ForceBase::setForces(const ParticleSet::ParticlePos& forces) { forces_ = forces; }

void ForceBase::setForces(Real val) { forces_ = val; }

void ForceBase::setForcesIonIon(const ParticleSet::ParticlePos& forces_ion_ion) { forces_ion_ion_ = forces_ion_ion; }

void ForceBase::initVarReduction(Real rcut, int m, int numFuncs)
{
  m_    = m;
  rcut_ = rcut;
  std::vector<Real> h(numFuncs);
  Matrix<Real> S(numFuncs, numFuncs);
  ck_.resize(numFuncs, 0.0);
  Real R2jp1 = rcut_ * rcut_;
  Real R2m   = 1.0;
  for (int i = 0; i < m_; i++)
    R2m *= rcut_;
  for (int j = 1; j <= numFuncs; j++)
  {
    h[j - 1] = R2jp1 / Real(j + 1);
    Real R2k = rcut_;
    for (int k = 1; k <= numFuncs; k++)
    {
      S(k - 1, j - 1) = R2m * R2k * R2jp1 / (Real)(m_ + k + j + 1);
      S(k - 1, j - 1) = std::pow(rcut_, (m_ + k + j + 1)) / (m_ + k + j + 1.0);
      R2k *= rcut_;
    }
    R2jp1 *= rcut_;
  }
  invert_matrix(S, false);
  for (int i = 0; i < numFuncs; i++)
  {
    for (int j = 0; j < numFuncs; j++)
      ck_[i] += S(i, j) * h[j];
  }
  FILE* fout = fopen("g_r.dat", "w");
  for (double r = 0.0; r < rcut_; r += 0.001)
    fprintf(fout, "%1.10f %1.10e\n", r, g(r));
  fclose(fout);
  app_log() << "Initialized variance reduction coefs.\n";
}

} // namespace qmcplusplus
