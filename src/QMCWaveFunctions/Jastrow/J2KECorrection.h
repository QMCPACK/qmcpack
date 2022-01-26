//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2021 QMCPACK developers.
//
// File developed by: Jeongnim Kim, jeongnim.kim@intel.com, Intel Corp.
//                    Amrita Mathuriya, amrita.mathuriya@intel.com, Intel Corp.
//                    Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//
// File created by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//////////////////////////////////////////////////////////////////////////////////////
// -*- C++ -*-
#ifndef QMCPLUSPLUS_J2KECORRECTION_H
#define QMCPLUSPLUS_J2KECORRECTION_H

#include <cmath>
#include <vector>
#include <ParticleSet.h>

namespace qmcplusplus
{
// helper class to activate KEcorr during optimizing Jastrow
template<typename RT, class FT>
class J2KECorrection
{
  size_t num_groups_;
  std::vector<size_t> num_elec_in_groups_;
  RT num_elecs_;
  RT vol;
  RT G0mag;
  const std::vector<FT*>& F_;
  bool SK_enabled;

public:
  J2KECorrection(const ParticleSet& targetPtcl, const std::vector<FT*>& F)
      : num_groups_(targetPtcl.groups()),
        num_elecs_(targetPtcl.getTotalNum()),
        vol(targetPtcl.getLattice().Volume),
        F_(F),
        SK_enabled(targetPtcl.hasSK())
  {
    // compute num_elec_in_groups_
    num_elec_in_groups_.reserve(3);
    for (int i = 0; i < num_groups_; i++)
      num_elec_in_groups_.push_back(targetPtcl.last(i) - targetPtcl.first(i));

    if (SK_enabled)
      G0mag = std::sqrt(targetPtcl.getSimulationCell().getKLists().ksq[0]);
  }

  RT computeKEcorr()
  {
    if (!SK_enabled)
      return 0;

    const int numPoints = 1000;
    RT uk               = 0.0;
    RT a                = 1.0;

    for (int i = 0; i < num_groups_; i++)
    {
      int Ni = num_elec_in_groups_[i];
      for (int j = 0; j < num_groups_; j++)
      {
        int Nj = num_elec_in_groups_[j];
        if (F_[i * num_groups_ + j])
        {
          FT& ufunc = *(F_[i * num_groups_ + j]);
          RT radius = ufunc.cutoff_radius;
          RT k      = G0mag;
          RT dr     = radius / (RT)(numPoints - 1);
          for (int ir = 0; ir < numPoints; ir++)
          {
            RT r = dr * (RT)ir;
            RT u = ufunc.evaluate(r);
            uk += 0.5 * 4.0 * M_PI * r * std::sin(k * r) / k * u * dr * (RT)Nj / (RT)(Ni + Nj);
          }
        }
      }
    }
    for (int iter = 0; iter < 20; iter++)
      a = uk / (4.0 * M_PI * (1.0 / (G0mag * G0mag) - 1.0 / (G0mag * G0mag + 1.0 / a)));
    return 4.0 * M_PI * a / (4.0 * vol) * num_elecs_;
  }
};
} // namespace qmcplusplus
#endif
