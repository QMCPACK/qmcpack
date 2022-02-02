//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2020 QMCPACK developers.
//
// File developed by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////


#include "EwaldHandler.h"

namespace qmcplusplus
{
void EwaldHandler::initBreakup(ParticleSet& ref)
{
  SuperCellEnum = ref.getLattice().SuperCellEnum;
  LR_rc         = ref.getLattice().LR_rc;
  LR_kc         = ref.getLattice().LR_kc;
  Sigma         = 3.5;
  //determine the sigma
  while (erfc(Sigma) / LR_rc > 1e-10)
  {
    Sigma += 0.1;
  }
  app_log() << "   EwaldHandler Sigma/LR_rc = " << Sigma;
  Sigma /= ref.getLattice().LR_rc;
  app_log() << "  Sigma=" << Sigma << std::endl;
  Volume     = ref.getLattice().Volume;
  PreFactors = 0.0;
  //See A.18
  if (SuperCellEnum == SUPERCELL_SLAB)
  {
    Area          = std::abs(ref.getLattice().R(0, 0) * ref.getLattice().R(1, 1) - ref.getLattice().R(0, 1) * ref.getLattice().R(1, 0));
    PreFactors[0] = 2.0 * M_PI / Area;                      //\f$ \frac{2\pi}{A}\f$
    PreFactors[1] = 2.0 * std::sqrt(M_PI) / (Sigma * Area); //\f$  \frac{2\pi}{A}\frac{1}{Sigma\pi}\f$
    PreFactors[2] = 1.0 / (2 * Sigma);                      //Used for the k-dependent term
  }
  fillFk(ref.getSimulationCell().getKLists());
}

EwaldHandler::EwaldHandler(const EwaldHandler& aLR, ParticleSet& ref)
    : LRHandlerBase(aLR), Sigma(aLR.Sigma), Volume(aLR.Volume), Area(aLR.Area), PreFactors(aLR.PreFactors)
{
  SuperCellEnum = aLR.SuperCellEnum;
  if (SuperCellEnum == SUPERCELL_SLAB)
    kMag = aLR.kMag;
}

void EwaldHandler::fillFk(const KContainer& KList)
{
  Fk.resize(KList.kpts_cart.size());
  const std::vector<int>& kshell(KList.kshell);
  MaxKshell = kshell.size() - 1;
  Fk_symm.resize(MaxKshell);
  kMag.resize(MaxKshell);
  if (SuperCellEnum == SUPERCELL_SLAB)
  {
    mRealType knorm           = M_PI / Area;
    mRealType ksum            = 0.0;
    mRealType oneovertwosigma = 1.0 / (2.0 * Sigma);
    for (int ks = 0, ki = 0; ks < Fk_symm.size(); ks++)
    {
      kMag[ks]     = std::sqrt(KList.ksq[ki]);
      mRealType uk = knorm / kMag[ks]; //pi/(A*k)
      Fk_symm[ks]  = uk;
      while (ki < KList.kshell[ks + 1] && ki < Fk.size())
        Fk[ki++] = uk;
      ksum += uk * erfc(kMag[ks] * oneovertwosigma) * (KList.kshell[ks + 1] - KList.kshell[ks]);
    }
    //correction for the self-energy returned by evaluateLR_r0
    PreFactors[3] = 2.0 * (std::sqrt(M_PI) / (Area * Sigma) - ksum);
  }
  else
  {
#if OHMMS_DIM == 2
    mRealType kgauss = 1.0 / (4 * Sigma * Sigma);
    mRealType knorm  = 2 * M_PI / Volume;
    for (int ks = 0, ki = 0; ks < Fk_symm.size(); ks++)
    {
      mRealType t2e = KList.ksq[ki] * kgauss;
      mRealType uk  = knorm * std::exp(-t2e) / KList.ksq[ki];
      Fk_symm[ks]   = uk;
      while (ki < KList.kshell[ks + 1] && ki < Fk.size())
        Fk[ki++] = uk;
    }
    PreFactors[3] = 0.0;
#elif OHMMS_DIM == 3
    mRealType kgauss = 1.0 / (4 * Sigma * Sigma);
    mRealType knorm  = 4 * M_PI / Volume;
    for (int ks = 0, ki = 0; ks < Fk_symm.size(); ks++)
    {
      mRealType t2e = KList.ksq[ki] * kgauss;
      mRealType uk  = knorm * std::exp(-t2e) / KList.ksq[ki];
      Fk_symm[ks]   = uk;
      while (ki < KList.kshell[ks + 1] && ki < Fk.size())
        Fk[ki++] = uk;
    }
    PreFactors[3] = 0.0;
#endif
  }
  app_log().flush();
}

EwaldHandler::mRealType EwaldHandler::evaluate_vlr_k(mRealType k) const
{
  mRealType uk = 0.0;
  if (SuperCellEnum == SUPERCELL_SLAB)
  {
    mRealType knorm = M_PI / Area;
    uk              = knorm / k; //pi/(A*k)
  }
  else
  {
#if OHMMS_DIM == 2
    mRealType kgauss = 1.0 / (4 * Sigma * Sigma);
    mRealType knorm  = 2 * M_PI / Volume;
    mRealType k2     = k * k;
    uk               = knorm * std::exp(-k2 * kgauss) / k2;
#elif OHMMS_DIM == 3
    mRealType kgauss = 1.0 / (4 * Sigma * Sigma);
    mRealType knorm  = 4 * M_PI / Volume;
    mRealType k2     = k * k;
    uk               = knorm * std::exp(-k2 * kgauss) / k2;
  }
#endif
    return uk;
  }

  EwaldHandler::mRealType EwaldHandler::evaluate_slab(pRealType z, const std::vector<int>& kshell,
                                                      const pComplexType* restrict eikr_i,
                                                      const pComplexType* restrict eikr_j) const
  {
    mRealType zp = z * Sigma;
    mRealType vk = -SlabFunc0(z, zp);
    //cout << "### SLAB " << z << " " << zp << std::endl;
    for (int ks = 0, ki = 0; ks < MaxKshell; ks++)
    {
      mRealType u = 0; //\sum Real (e^ikr_i e^(-ikr_j))
      for (; ki < kshell[ks + 1]; ki++, eikr_i++, eikr_j++)
        u += ((*eikr_i).real() * (*eikr_j).real() + (*eikr_i).imag() * (*eikr_j).imag());
      vk += u * Fk_symm[ks] * SlabFuncK(ks, z, zp);
    }
    return vk;
  }
} // namespace qmcplusplus
