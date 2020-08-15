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

#ifndef QMCPLUSPLUS_AFQMC_WALKERUPDATES_HPP
#define QMCPLUSPLUS_AFQMC_WALKERUPDATES_HPP

#include <algorithm>

#include "Configuration.h"
#include "AFQMC/config.h"

namespace qmcplusplus
{
namespace afqmc
{
template<class Wlk>
void free_projection_walker_update(Wlk&& w,
                                   RealType dt,
                                   ComplexType overlap,
                                   ComplexType MFfactor,
                                   RealType Eshift,
                                   ComplexType hybrid_weight)
{
  ComplexType old_ovlp = *w.overlap();
  ComplexType old_eloc = *w.pseudo_energy();
  ComplexType eloc;
  RealType scale            = 1.0;
  ComplexType ratioOverlaps = ComplexType(1.0, 0.0);
  eloc                      = MFfactor / dt;
  ComplexType factor        = std::exp(-dt * (0.5 * (eloc + old_eloc) - Eshift));
  *w.weight() *= std::abs(factor);
  *w.phase() *= factor / std::abs(factor);
  *w.pseudo_energy() = eloc;
  *w.overlap()       = overlap;
}

template<class Wlk>
void hybrid_walker_update(Wlk&& w,
                          RealType dt,
                          bool apply_constrain,
                          bool imp_sampl,
                          RealType Eshift,
                          ComplexType overlap,
                          ComplexType MFfactor,
                          ComplexType hybrid_weight)
{
  ComplexType old_ovlp = *w.overlap();
  ComplexType old_eloc = *w.pseudo_energy();
  ComplexType eloc;
  RealType scale            = 1.0;
  ComplexType ratioOverlaps = ComplexType(1.0, 0.0);

  if (imp_sampl)
    ratioOverlaps = overlap / old_ovlp;
  if (!std::isfinite(ratioOverlaps.real()) && apply_constrain && imp_sampl)
  {
    scale = 0.0;
    eloc  = old_eloc;
  }
  else
  {
    scale = (apply_constrain ? (std::max(0.0, std::cos(std::arg(ratioOverlaps) - MFfactor.imag()))) : 1.0);
    if (imp_sampl)
      eloc = (MFfactor - hybrid_weight - std::log(ratioOverlaps)) / dt;
    else
      eloc = MFfactor / dt;
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

#if defined(__AFQMC_DEBUG_VERBOSE__)
  std::cout << " update: "
            << "    eloc:          " << eloc << "\n"
            << "    eloc_:         " << eloc_ << "\n"
            << "    ov:            " << overlap << "\n"
            << "    old_ov:        " << old_ovlp << "\n"
            << "    ratio:         " << ratioOverlaps << "\n"
            << "    MFfactor:      " << MFfactor << "\n"
            << "    hybrid_weight: " << hybrid_weight << "\n"
            << "    scale:         " << scale << "\n"
            << std::endl;
#endif

  *w.weight() *= ComplexType(scale * std::exp(-dt * (0.5 * (eloc.real() + old_eloc.real()) - Eshift)), 0.0);
  if (w.NumBackProp() > 0 && std::abs(scale) > 1e-16)
  {
    *w.BPWeightFactor() *= std::exp(-ComplexType(0.0, dt) * (0.5 * (eloc.imag() + old_eloc.imag()))) / scale;
  }

  *w.pseudo_energy() = eloc;
  *w.overlap()       = overlap;
}

template<class Wlk, class TVec>
void local_energy_walker_update(Wlk&& w,
                                RealType dt,
                                bool apply_constrain,
                                RealType Eshift,
                                ComplexType overlap,
                                TVec&& energies,
                                ComplexType MFfactor,
                                ComplexType hybrid_weight)
{
  ComplexType old_ovlp      = *w.overlap();
  ComplexType old_eloc      = *w.pseudo_energy();
  ComplexType eloc          = energies[0] + energies[1] + energies[2];
  RealType scale            = 1.0;
  ComplexType ratioOverlaps = overlap / old_ovlp;

  if (!std::isfinite((ratioOverlaps * MFfactor).real()) && apply_constrain)
  {
    scale = 0.0;
    eloc  = old_eloc;
  }
  else
    scale = (apply_constrain ? (std::max(0.0, std::cos(std::arg(ratioOverlaps) - MFfactor.imag()))) : 1.0);

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

  *w.weight() *= ComplexType(scale * std::exp(-dt * (0.5 * (eloc.real() + old_eloc.real()) - Eshift)), 0.0);
  if (w.NumBackProp() > 0 && std::abs(scale) > 1e-16)
  {
    *w.BPWeightFactor() *= std::exp(-ComplexType(0.0, dt) * (0.5 * (eloc.imag() + old_eloc.imag()))) / scale;
  }
  *w.pseudo_energy() = eloc;
  *w.E1()            = energies[0];
  *w.EXX()           = energies[1];
  *w.EJ()            = energies[2];
  *w.overlap()       = overlap;
}

} // namespace afqmc

} // namespace qmcplusplus

#endif
