//////////////////////////////////////////////////////////////////
// (c) Copyright 2003-  by Jeongnim Kim and Kenneth P Esler
//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////
//   National Center for Supercomputing Applications &
//   Materials Computation Center
//   University of Illinois, Urbana-Champaign
//   Urbana, IL 61801
//   e-mail: jnkim@ncsa.uiuc.edu
//
// Supported by
//   National Center for Supercomputing Applications, UIUC
//   Materials Computation Center, UIUC
//////////////////////////////////////////////////////////////////
// -*- C++ -*-
#include "QMCWaveFunctions/OrbitalBase.h"
#include "QMCWaveFunctions/Fermion/SlaterDet.h"

namespace qmcplusplus
{

// void
// ratio (MCWalkerConfiguration &W, int iat, vector<PosType> &new_pos,
// 	   vector<ValueType> &psi_ratios)
// {
//   Dets[DetID[iat]]->ratio(W, iat, new_pos, psi_ratios);
// }

// void
// ratio (MCWalkerConfiguration &W, int iat, vector<PosType> &new_pos,
// 	   vector<ValueType> &psi_ratios,	vector<GradType>  &grad)
// {
//   Dets[DetID[iat]]->ratio(W, iat, new_pos, psi_ratios, grad);
// }

void
SlaterDet::ratio (vector<Walker_t*> &walkers,    vector<int> &iatList,
                  vector<PosType> &rNew, vector<ValueType> &psi_ratios,
                  vector<GradType>  &grad, vector<ValueType> &lapl)
{
  // Sort walkers by determinant number
  vector<vector<Walker_t*> > sorted_walkers(Dets.size());
  vector<vector<int> >       sorted_iatList(Dets.size());
  vector<vector<PosType> >   sorted_rNew(Dets.size());
  vector<vector<ValueType> > ratio_det(Dets.size()), lapl_det(Dets.size());
  vector<vector<GradType> >  grad_det(Dets.size());
  for (int iw=0; iw<walkers.size(); iw++)
  {
    int det = DetID[iatList[iw]];
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
  vector<int> index(Dets.size());
  for (int iw=0; iw<walkers.size(); iw++)
  {
    int det = DetID[iatList[iw]];
    int i = index[det]++;
    psi_ratios[iw] = ratio_det[det][i];
    grad[iw] = grad_det[det][i];
    lapl[iw] = lapl_det[det][i];
  }
}

void SlaterDet::update (const vector<Walker_t*> &walkers,
                        const vector<int> &iatList)
{
  // Sort walkers by determinant number
  vector<vector<Walker_t*> > sorted_walkers(Dets.size());
  vector<vector<int> >       sorted_iatList(Dets.size());
  for (int iw=0; iw<walkers.size(); iw++)
  {
    int det = DetID[iatList[iw]];
    sorted_walkers[det].push_back(walkers[iw]);
    sorted_iatList[det].push_back(iatList[iw]);
  }
  // Call each DiracDeterminant with the appropriate walkers
  for (int idet=0; idet<Dets.size(); idet++)
    if (sorted_walkers[idet].size())
      Dets[idet]->update(sorted_walkers[idet], sorted_iatList[idet]);
}

}
/***************************************************************************
 * $RCSfile$   $Author: jmcminis $
 * $Revision: 4646 $   $Date: 2010-02-23 00:25:45 -0600 (Tue, 23 Feb 2010) $
 * $Id: SlaterDet.h 4646 2010-02-23 06:25:45Z jmcminis $
 ***************************************************************************/
