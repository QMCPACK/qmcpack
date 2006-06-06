#ifndef QMCPLUSPLUS_DISPERSIONRELATION_H
#define QMCPLUSPLUS_DISPERSIONRELATION_H

#include "Configuration.h"
#include "LongRange/KContainer.h"
#include "Particle/ParticleSet.h"
#include "OhmmsPETE/OhmmsVector.h"
#include "OhmmsPETE/TinyVector.h"

namespace qmcplusplus {
  typedef QMCTraits::RealType RealType;
  typedef QMCTraits::PosType  PosType;
  typedef QMCTraits::IndexType IndexType;
  

  void gen3DLattice(RealType cut, ParticleSet& ptclSet, 
                    Vector<PosType>& LocLattice, Vector<IndexType>& DimSizes) {
    KContainer localContainer(ptclSet.Lattice);
    localContainer.UpdateKLists(cut, false);
    LocLattice.resize(localContainer.kpts_cart.size());
    for (IndexType i = 0; i < LocLattice.size(); i++) {
      LocLattice[i] = localContainer.kpts_cart[i];
    }
    DimSizes.resize(3);
    for (IndexType i = 0; i < 3; i++) {
      DimSizes[i] = localContainer.mmax[i] * 2;
    }
  }
  
  void gen1DLattice(RealType cut, const ParticleSet& ptclSet, 
                    Vector<PosType>& LocLattice, Vector<IndexType>& DimSizes) {
    RealType GVectorMomentum = std::sqrt(dot(ptclSet.Lattice.k_cart(PosType(1,0,0)), ptclSet.Lattice.k_cart(PosType(1,0,0))));
    RealType CutoffMomentum = std::sqrt(2 * cut);
    DimSizes.resize(1);
    DimSizes[0] = int(std::floor(CutoffMomentum / GVectorMomentum)) * 2;
    LocLattice.resize(DimSizes[0]);
    for (IndexType i = 0; i < DimSizes[0]; i++) {
      IndexType multiple = i;
      if (multiple > DimSizes[0]/2) multiple -= DimSizes[0];
      LocLattice[i] = multiple * ptclSet.Lattice.k_cart(PosType(1,0,0));
    }    
  }
  
  void genFreeParticleDispersion(const Vector<PosType>& LocLattice, Vector<RealType>& Disp) {
    Disp.resize(LocLattice.size());
    for (IndexType i = 0; i < LocLattice.size(); i++) {
      Disp[i] = dot(LocLattice[i], LocLattice[i]) / 2.0;
    }
  }

  void genSimpleModelDispersion(const Vector<PosType>& LocLattice, 
                                Vector<RealType>& Disp, RealType GapSize,
                                RealType FermiMomentum) {
    Disp.resize(LocLattice.size());
    RealType FermiEnergy = FermiMomentum * FermiMomentum / 2.0;
    for (IndexType i = 0; i < LocLattice.size(); i++) {
      RealType temp = dot(LocLattice[i], LocLattice[i]) / 2.0;
      if (temp > FermiEnergy) {
        Disp[i] = temp + GapSize;
      } else {
        Disp[i] = temp;
      }
    }
  }
    
  void genPennModelDispersion(const Vector<PosType>& LocLattice, 
                                Vector<RealType>& Disp, RealType GapSize,
                                RealType FermiMomentum) {
    Disp.resize(LocLattice.size());
    RealType FermiEnergy = FermiMomentum * FermiMomentum / 2.0;
    RealType GapSq = GapSize * GapSize;
    for (IndexType i = 0; i < LocLattice.size(); i++) {
      RealType currMom = std::sqrt(dot(LocLattice[i], LocLattice[i]));
      RealType kEnergy = currMom * currMom / 2.0;
      RealType kPrimeEnergy = (currMom - 2 * FermiMomentum) * (currMom - 2 * FermiMomentum) / 2.0;
      if (currMom > FermiMomentum) {
        Disp[i] = 0.5 * (kEnergy + kPrimeEnergy +
                         std::sqrt(std::pow(kEnergy - kPrimeEnergy, 2) + GapSq));
      } else {
        Disp[i] = 0.5 * (kEnergy + kPrimeEnergy - 
                         std::sqrt(std::pow(kEnergy - kPrimeEnergy, 2) + GapSq));
      }
    }
  }
  
}

#endif
    
