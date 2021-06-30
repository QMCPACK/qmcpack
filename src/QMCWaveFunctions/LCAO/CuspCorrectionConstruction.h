//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Mark Dewing, mdewing@anl.gov, Argonne National Laboratory
//                    Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//
// File created by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//////////////////////////////////////////////////////////////////////////////////////


#ifndef QMCPLUSPLUS_CUSP_CORRECTION_CONSTRUCTOR_H
#define QMCPLUSPLUS_CUSP_CORRECTION_CONSTRUCTOR_H

#include "LCAOrbitalSet.h"
#include "LCAOrbitalSetWithCorrection.h"
#include "CuspCorrection.h"

namespace qmcplusplus
{
// Modifies orbital set lcwc
void applyCuspCorrection(const Matrix<CuspCorrectionParameters>& info,
                         int num_centers,
                         int orbital_set_size,
                         ParticleSet& targetPtcl,
                         ParticleSet& sourcePtcl,
                         LCAOrbitalSetWithCorrection& lcwc,
                         const std::string& id);

void saveCusp(int orbital_set_size, int num_centers, Matrix<CuspCorrectionParameters>& info, const std::string& id);

void generateCuspInfo(int orbital_set_size,
                      int num_centers,
                      Matrix<CuspCorrectionParameters>& info,
                      const ParticleSet& targetPtcl,
                      const ParticleSet& sourcePtcl,
                      const LCAOrbitalSetWithCorrection& lcwc,
                      const std::string& id,
                      Communicate& Comm);
} // namespace qmcplusplus

#endif
