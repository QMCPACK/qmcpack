//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////
    
    
#include "DispersionRelation.h"
#include "LongRange/KContainer.h"

namespace qmcplusplus
{
typedef QMCTraits::RealType   RealType;
typedef QMCTraits::IndexType  IndexType;
typedef QMCTraits::PosType    PosType;

void gen3DLattice(RealType cut, ParticleSet& ptclSet,
                  Vector<PosType>& LocLattice, Vector<RealType>& Degeneracies,
                  Vector<IndexType>& DimSizes)
{
  KContainer localContainer(ptclSet.Lattice);
  localContainer.UpdateKLists(cut, false);
  LocLattice.resize(localContainer.kpts_cart.size());
  Degeneracies.resize(localContainer.kpts_cart.size());
  for (IndexType i = 0; i < LocLattice.size(); i++)
  {
    LocLattice[i] = localContainer.kpts_cart[i];
    Degeneracies[i] = 1.0;
  }
  DimSizes.resize(3);
  for (IndexType i = 0; i < 3; i++)
  {
    DimSizes[i] = localContainer.mmax[i] * 2;
  }
}

void gen1DLattice(RealType cut, ParticleSet& ptclSet,
                  Vector<PosType>& LocLattice, Vector<RealType>& Degeneracies,
                  Vector<IndexType>& DimSizes)
{
  RealType GVectorMomentum = std::sqrt(dot(ptclSet.Lattice.k_cart(PosType(1,0,0)), ptclSet.Lattice.k_cart(PosType(1,0,0))));
  RealType CutoffMomentum = std::sqrt(2 * cut);
  DimSizes.resize(1);
  DimSizes[0] = int(std::floor(CutoffMomentum / GVectorMomentum)) * 2;
  LocLattice.resize(DimSizes[0]);
  Degeneracies.resize(DimSizes[0]);
  for (IndexType i = 0; i < DimSizes[0]; i++)
  {
    IndexType multiple = i;
    Degeneracies[i] = 1.0;
    if (multiple > DimSizes[0]/2)
      multiple -= DimSizes[0];
    LocLattice[i] = multiple * ptclSet.Lattice.k_cart(PosType(1,0,0));
  }
}

void genDegenLattice(RealType cut, ParticleSet& ptclSet,
                     Vector<PosType>& LocLattice, Vector<RealType>& Degeneracies,
                     Vector<IndexType>& DimSizes)
{
  KContainer localContainer(ptclSet.Lattice);
  localContainer.UpdateKLists(cut, true);
  LocLattice.resize(localContainer.kpts_sorted.size()+1);
  Degeneracies.resize(localContainer.kpts_sorted.size()+1);
  LocLattice[0] = PosType(0,0,0);
  Degeneracies[0] = 1.0;
  std::cout << "In genDegenLattice, about to start loop over kpts_sorted." << std::endl;
  std::map<int, std::vector<int>*>::iterator it(localContainer.kpts_sorted.begin());
  for (IndexType i = 1; i < LocLattice.size(); i++)
  {
    std::vector<int>& vecref(*(it->second));
    Degeneracies[i] = vecref.size();
    LocLattice[i] = localContainer.kpts_cart[vecref[0]];
    ++it;
  }
  for (IndexType i = 0; i < LocLattice.size(); i++)
  {
    std::cout << Degeneracies[i] << "      " << LocLattice[i] << std::endl;
  }
}

void genFreeParticleDispersion(const Vector<PosType>& LocLattice, Vector<RealType>& Disp)
{
  for (IndexType i = 0; i < LocLattice.size(); i++)
  {
    Disp[i] *= dot(LocLattice[i], LocLattice[i]) / 2.0;
  }
}

void genDebugDispersion(const Vector<PosType>& LocLattice, Vector<RealType>& Disp)
{
  for (IndexType i = 0; i < LocLattice.size(); i++)
  {
    Disp[i] = 1.0;
  }
  Disp[0] = 1.0;
}

void genSimpleModelDispersion(const Vector<PosType>& LocLattice,
                              Vector<RealType>& Disp, RealType GapSize,
                              RealType FermiMomentum)
{
  RealType FermiEnergy = FermiMomentum * FermiMomentum / 2.0;
  for (IndexType i = 0; i < LocLattice.size(); i++)
  {
    RealType temp = dot(LocLattice[i], LocLattice[i]) / 2.0;
    if (temp > FermiEnergy)
    {
      Disp[i] *= temp + GapSize;
    }
    else
    {
      Disp[i] *= temp;
    }
  }
}

void genPennModelDispersion(const Vector<PosType>& LocLattice,
                            Vector<RealType>& Disp, RealType GapSize,
                            RealType FermiMomentum)
{
  RealType FermiEnergy = FermiMomentum * FermiMomentum / 2.0;
  RealType GapSq = GapSize * GapSize;
  for (IndexType i = 0; i < LocLattice.size(); i++)
  {
    RealType currMom = std::sqrt(dot(LocLattice[i], LocLattice[i]));
    RealType kEnergy = currMom * currMom / 2.0;
    RealType kPrimeEnergy = (currMom - 2 * FermiMomentum) * (currMom - 2 * FermiMomentum) / 2.0;
    if (currMom > FermiMomentum)
    {
      Disp[i] *= 0.5 * (kEnergy + kPrimeEnergy +
                        std::sqrt(std::pow(kEnergy - kPrimeEnergy, 2) + GapSq));
    }
    else
    {
      Disp[i] *= 0.5 * (kEnergy + kPrimeEnergy -
                        std::sqrt(std::pow(kEnergy - kPrimeEnergy, 2) + GapSq));
    }
  }
}
}
