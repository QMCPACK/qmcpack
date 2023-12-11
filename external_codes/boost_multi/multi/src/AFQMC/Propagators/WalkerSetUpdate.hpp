//////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source
// License.  See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by:
// Miguel A. Morales, moralessilva2@llnl.gov
//    Lawrence Livermore National Laboratory
//
// File created by:
// Miguel A. Morales, moralessilva2@llnl.gov
//    Lawrence Livermore National Laboratory
////////////////////////////////////////////////////////////////////////////////

#ifndef QMCPLUSPLUS_AFQMC_WALKERSETUPDATES_HPP
#define QMCPLUSPLUS_AFQMC_WALKERSETUPDATES_HPP

#include <algorithm>

#include "Configuration.h"
#include "AFQMC/config.h"
#include "AFQMC/Numerics/ma_blas.hpp"
#include "AFQMC/Walkers/WalkerConfig.hpp"

namespace qmcplusplus
{
namespace afqmc
{
template<class Wlk, class OMat, class Mat, class WMat>
void free_projection_walker_update(Wlk&& w,
                                   RealType dt,
                                   OMat&& overlap,
                                   Mat&& MFfactor,
                                   RealType Eshift,
                                   Mat&& hybrid_weight,
                                   WMat& work)
{
  int nwalk = w.size();
  // constexpr if can be used to avoid the memory copy, by comparing the pointer types
  // between WMat and Mat/OMat
  if (std::get<0>(work.sizes()) < 7 || std::get<1>(work.sizes()) < nwalk)
    work.reextent({7, nwalk});

  w.getProperty(WEIGHT, work[0]);
  w.getProperty(PHASE, work[1]);
  w.getProperty(PSEUDO_ELOC_, work[2]);
  w.getProperty(OVLP, work[3]);
  using std::copy_n;
  copy_n(overlap.origin(), nwalk, work[4].origin());
  copy_n(MFfactor.origin(), nwalk, work[5].origin());
  copy_n(hybrid_weight.origin(), nwalk, work[6].origin());

  for (int i = 0; i < nwalk; i++)
  {
    ComplexType old_ovlp = work[3][i];
    ComplexType old_eloc = work[2][i];
    ComplexType eloc;
    ComplexType ratioOverlaps = ComplexType(1.0, 0.0);
    eloc                      = work[5][i] / dt;
    ComplexType factor        = std::exp(-dt * (0.5 * (eloc + old_eloc) - Eshift));
    work[0][i] *= std::abs(factor);
    work[1][i] *= factor / std::abs(factor);
    work[2][i] = eloc;
    work[3][i] = work[4][i];
  }

  w.setProperty(WEIGHT, work[0]);
  w.setProperty(PHASE, work[1]);
  w.setProperty(PSEUDO_ELOC_, work[2]);
  w.setProperty(OVLP, work[3]);
}

template<class Wlk, class OMat, class Mat, class WMat>
void hybrid_walker_update(Wlk&& w,
                          RealType dt,
                          bool apply_constrain,
                          bool imp_sampl,
                          RealType Eshift,
                          OMat&& overlap,
                          Mat&& MFfactor,
                          Mat&& hybrid_weight,
                          WMat& work)
{
  int nwalk = w.size();
  // constexpr if can be used to avoid the memory copy, by comparing the pointer types
  // between WMat and Mat/OMat
  if (std::get<0>(work.sizes()) < 7 || std::get<1>(work.sizes()) < nwalk)
    work.reextent({7, nwalk});

  bool BackProp = (w.getBPPos() >= 0 && w.getBPPos() < w.NumBackProp());

  w.getProperty(WEIGHT, work[0]);
  w.getProperty(PSEUDO_ELOC_, work[1]);
  w.getProperty(OVLP, work[2]);
  using std::copy_n;
  copy_n(overlap.origin(), nwalk, work[4].origin());
  copy_n(MFfactor.origin(), nwalk, work[5].origin());
  copy_n(hybrid_weight.origin(), nwalk, work[6].origin());

  for (int i = 0; i < nwalk; i++)
  {
    ComplexType old_ovlp = work[2][i];
    ComplexType old_eloc = work[1][i];
    ComplexType eloc;
    RealType scale            = 1.0;
    ComplexType ratioOverlaps = ComplexType(1.0, 0.0);

    if (imp_sampl)
      ratioOverlaps = work[4][i] / old_ovlp;
    if (!std::isfinite(ratioOverlaps.real()) && apply_constrain && imp_sampl)
    {
      scale = 0.0;
      eloc  = old_eloc;
    }
    else
    {
      scale = (apply_constrain ? (std::max(0.0, std::cos(std::arg(ratioOverlaps) - work[5][i].imag()))) : 1.0);
      if (imp_sampl)
        eloc = (work[5][i] - work[6][i] - std::log(ratioOverlaps)) / dt;
      else
        eloc = work[5][i] / dt;
    }
    ComplexType eloc_ = eloc;

    if ((!std::isfinite(eloc.real())) || (std::abs(eloc.real()) < std::numeric_limits<RealType>::min()))
    {
      scale = 0.0;
      eloc  = old_eloc;
    }
    else
    {
      eloc = ComplexType(std::max(std::min(eloc.real(), Eshift + std::sqrt(2.0 / dt)), Eshift - std::sqrt(2.0 / dt)),
                         eloc.imag());
    }
//#define __AFQMC_DEBUG_VERBOSE__
#if defined(__AFQMC_DEBUG_VERBOSE__)
    std::cout << " update: "
              << "    eloc:          " << eloc << "\n"
              << "    eloc_:         " << eloc_ << "\n"
              << "    ov:            " << work[4][i] << "\n"
              << "    old_ov:        " << old_ovlp << "\n"
              << "    ratio:         " << ratioOverlaps << "\n"
              << "    MFfactor:      " << work[5][i] << "\n"
              << "    hybrid_weight: " << work[6][i] << "\n"
              << "    scale:         " << scale << "\n"
              << std::endl;
#endif

    work[0][i] *= ComplexType(scale * std::exp(-dt * (0.5 * (eloc.real() + old_eloc.real()) - Eshift)), 0.0);
    work[1][i] = eloc;
    work[2][i] = work[4][i];
    work[3][i] = std::exp(-ComplexType(0.0, dt) * (0.5 * (eloc.imag() + old_eloc.imag()))) / scale;
  }
  w.setProperty(WEIGHT, work[0]);
  w.setProperty(PSEUDO_ELOC_, work[1]);
  w.setProperty(OVLP, work[2]);
  if (BackProp)
  {
    auto&& WFac((*w.getWeightFactors())[w.getHistoryPos()]);
    using std::copy_n;
    copy_n(work[3].origin(), nwalk, WFac.origin());
    auto&& WHis((*w.getWeightHistory())[w.getHistoryPos()]);
    copy_n(work[0].origin(), nwalk, WHis.origin());
  }
}

template<class Wlk, class EMat, class OMat, class Mat, class WMat>
void local_energy_walker_update(Wlk&& w,
                                RealType dt,
                                bool apply_constrain,
                                RealType Eshift,
                                OMat&& overlap,
                                EMat&& energies,
                                Mat&& MFfactor,
                                Mat&& hybrid_weight,
                                WMat& work)
{
  int nwalk = w.size();
  // constexpr if can be used to avoid the memory copy, by comparing the pointer types
  // between WMat and Mat/OMat
  if (std::get<0>(work.sizes()) < 12 || std::get<1>(work.sizes()) < nwalk)
    work.reextent({12, nwalk});

  bool BackProp = (w.getBPPos() >= 0 && w.getBPPos() < w.NumBackProp());

  w.getProperty(WEIGHT, work[0]);
  w.getProperty(PSEUDO_ELOC_, work[1]);
  w.getProperty(OVLP, work[2]);
  w.getProperty(E1_, work[3]);
  w.getProperty(EXX_, work[4]);
  w.getProperty(EJ_, work[5]);
  using std::copy_n;
  copy_n(overlap.origin(), nwalk, work[7].origin());
  copy_n(MFfactor.origin(), nwalk, work[8].origin());
  ma::copy(energies({0, nwalk}, 0), work[9]);
  ma::copy(energies({0, nwalk}, 1), work[10]);
  ma::copy(energies({0, nwalk}, 2), work[11]);

  for (int i = 0; i < nwalk; i++)
  {
    ComplexType old_ovlp      = work[2][i];
    ComplexType old_eloc      = work[1][i];
    ComplexType eloc          = work[9][i] + work[10][i] + work[11][i];
    RealType scale            = 1.0;
    ComplexType ratioOverlaps = work[7][i] / old_ovlp;

    if (!std::isfinite((ratioOverlaps * work[8][i]).real()) && apply_constrain)
    {
      scale = 0.0;
      eloc  = old_eloc;
    }
    else
      scale = (apply_constrain ? (std::max(0.0, std::cos(std::arg(ratioOverlaps) - work[8][i].imag()))) : 1.0);

    if ((!std::isfinite(eloc.real())) || (std::abs(eloc.real()) < std::numeric_limits<RealType>::min()))
    {
      scale = 0.0;
      eloc  = old_eloc;
    }
    else
    {
      eloc = ComplexType(std::max(std::min(eloc.real(), Eshift + std::sqrt(2.0 / dt)), Eshift - std::sqrt(2.0 / dt)),
                         eloc.imag());
    }

    work[0][i] *= ComplexType(scale * std::exp(-dt * (0.5 * (eloc.real() + old_eloc.real()) - Eshift)), 0.0);
    work[6][i] = std::exp(-ComplexType(0.0, dt) * (0.5 * (eloc.imag() + old_eloc.imag()))) / scale;
    work[1][i] = eloc;
    work[2][i] = work[7][i];
    work[3][i] = work[9][i];
    work[4][i] = work[10][i];
    work[5][i] = work[11][i];
  }

  w.setProperty(WEIGHT, work[0]);
  w.setProperty(PSEUDO_ELOC_, work[1]);
  w.setProperty(OVLP, work[2]);
  w.setProperty(E1_, work[3]);
  w.setProperty(EXX_, work[4]);
  w.setProperty(EJ_, work[5]);
  if (BackProp)
  {
    auto&& WFac((*w.getWeightFactors())[w.getHistoryPos()]);
    using std::copy_n;
    copy_n(work[6].origin(), nwalk, WFac.origin());
    auto&& WHis((*w.getWeightHistory())[w.getHistoryPos()]);
    copy_n(work[0].origin(), nwalk, WHis.origin());
  }
}

} // namespace afqmc

} // namespace qmcplusplus

#endif
