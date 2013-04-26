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
