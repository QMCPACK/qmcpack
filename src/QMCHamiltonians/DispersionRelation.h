//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////
    
    
#ifndef QMCPLUSPLUS_DISPERSIONRELATION_H
#define QMCPLUSPLUS_DISPERSIONRELATION_H

#include "Configuration.h"
#include "Particle/ParticleSet.h"
#include "OhmmsPETE/OhmmsVector.h"

namespace qmcplusplus
{

void gen3DLattice(QMCTraits::RealType cut, ParticleSet& ptclSet,
                  Vector<QMCTraits::PosType>& LocLattice,
                  Vector<QMCTraits::RealType>& Degeneracies,
                  Vector<QMCTraits::IndexType>& DimSizes);

void gen1DLattice(QMCTraits::RealType cut, ParticleSet& ptclSet,
                  Vector<QMCTraits::PosType>& LocLattice,
                  Vector<QMCTraits::RealType>& Degeneracies,
                  Vector<QMCTraits::IndexType>& DimSizes);

void genDegenLattice(QMCTraits::RealType cut, ParticleSet& ptclSet,
                     Vector<QMCTraits::PosType>& LocLattice,
                     Vector<QMCTraits::RealType>& Degeneracies,
                     Vector<QMCTraits::IndexType>& DimSizes);

void genFreeParticleDispersion(const Vector<QMCTraits::PosType>& LocLattice,
                               Vector<QMCTraits::RealType>& Disp);

void genDebugDispersion(const Vector<QMCTraits::PosType>& LocLattice,
                        Vector<QMCTraits::RealType>& Disp);

void genSimpleModelDispersion(const Vector<QMCTraits::PosType>& LocLattice,
                              Vector<QMCTraits::RealType>& Disp,
                              QMCTraits::RealType GapSize,
                              QMCTraits::RealType FermiMomentum);

void genPennModelDispersion(const Vector<QMCTraits::PosType>& LocLattice,
                            Vector<QMCTraits::RealType>& Disp,
                            QMCTraits::RealType GapSize,
                            QMCTraits::RealType FermiMomentum);

}

#endif
