//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Ken Esler, kpesler@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////
    
    
#include "QMCWaveFunctions/OrbitalBase.h"
#include "QMCWaveFunctions/Fermion/SlaterDet.h"

namespace qmcplusplus
{

// void
// ratio (MCWalkerConfiguration &W, int iat, std::vector<PosType> &new_pos,
// 	   std::vector<ValueType> &psi_ratios)
// {
//   Dets[DetID[iat]]->ratio(W, iat, new_pos, psi_ratios);
// }

// void
// ratio (MCWalkerConfiguration &W, int iat, std::vector<PosType> &new_pos,
// 	   std::vector<ValueType> &psi_ratios,	vector<GradType>  &grad)
// {
//   Dets[DetID[iat]]->ratio(W, iat, new_pos, psi_ratios, grad);
// }

void
SlaterDet::ratio (std::vector<Walker_t*> &walkers,    std::vector<int> &iatList,
                  std::vector<PosType> &rNew, std::vector<ValueType> &psi_ratios,
                  std::vector<GradType>  &grad, std::vector<ValueType> &lapl)
{
  // Sort walkers by determinant number
  std::vector<std::vector<Walker_t*> > sorted_walkers(Dets.size());
  std::vector<std::vector<int> >       sorted_iatList(Dets.size());
  std::vector<std::vector<PosType> >   sorted_rNew(Dets.size());
  std::vector<std::vector<ValueType> > ratio_det(Dets.size()), lapl_det(Dets.size());
  std::vector<std::vector<GradType> >  grad_det(Dets.size());
  for (int iw=0; iw<walkers.size(); iw++)
  {
    int det = getDetID(iatList[iw]);
    sorted_walkers[det].push_back(walkers[iw]);
    sorted_iatList[det].push_back(iatList[iw]);
    sorted_rNew[det].push_back(rNew[iw]);
  }
  // Call each DiracDeterminant with the appropriate walkers
  for (int idet=0; idet<Dets.size(); idet++)
  {
    ratio_det[idet].resize(sorted_walkers[idet].size());
    grad_det[idet].resize(sorted_walkers[idet].size());
    lapl_det[idet].resize(sorted_walkers[idet].size());
    if (sorted_walkers[idet].size())
      Dets[idet]->ratio(sorted_walkers[idet], sorted_iatList[idet], sorted_rNew[idet],
                        ratio_det[idet], grad_det[idet], lapl_det[idet]);
  }
  // Copy ratios back into output
  std::vector<int> index(Dets.size());
  for (int iw=0; iw<walkers.size(); iw++)
  {
    int det = getDetID(iatList[iw]);
    int i = index[det]++;
    psi_ratios[iw] = ratio_det[det][i];
    grad[iw] = grad_det[det][i];
    lapl[iw] = lapl_det[det][i];
  }
}

void SlaterDet::update (const std::vector<Walker_t*> &walkers,
                        const std::vector<int> &iatList)
{
  // Sort walkers by determinant number
  std::vector<std::vector<Walker_t*> > sorted_walkers(Dets.size());
  std::vector<std::vector<int> >       sorted_iatList(Dets.size());
  for (int iw=0; iw<walkers.size(); iw++)
  {
    int det = getDetID(iatList[iw]);
    sorted_walkers[det].push_back(walkers[iw]);
    sorted_iatList[det].push_back(iatList[iw]);
  }
  // Call each DiracDeterminant with the appropriate walkers
  for (int idet=0; idet<Dets.size(); idet++)
    if (sorted_walkers[idet].size())
      Dets[idet]->update(sorted_walkers[idet], sorted_iatList[idet]);
}

}
