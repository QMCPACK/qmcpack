//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Bryan Clark, bclark@Princeton.edu, Princeton University
//                    Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Ken Esler, kpesler@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Bryan Clark, bclark@Princeton.edu, Princeton University
//////////////////////////////////////////////////////////////////////////////////////
    
    
#ifndef QMCPLUSPLUS_DENSITY_HAMILTONIAN_H
#define QMCPLUSPLUS_DENSITY_HAMILTONIAN_H
#include <QMCHamiltonians/QMCHamiltonianBase.h>
#include <OhmmsPETE/OhmmsArray.h>
#include "LongRange/LRCoulombSingleton.h"
namespace qmcplusplus
{

typedef LRCoulombSingleton::LRHandlerType LRHandlerType;
typedef LRCoulombSingleton::GridType       GridType;
typedef LRCoulombSingleton::RadFunctorType RadFunctorType;

class DensityEstimator: public QMCHamiltonianBase
{
public:

  DensityEstimator(ParticleSet& elns);
  int potentialIndex;
  void resetTargetParticleSet(ParticleSet& P);

  ///For Potential
  RealType evalSR(ParticleSet& P,int ipart);
  RealType evalLR(ParticleSet& P,int iat);
  void InitPotential(ParticleSet &P);
  std::vector<RealType> Zat,Zspec;
  RadFunctorType* rVs;
  int NumSpecies;
  int NumCenters;
  LRHandlerType* AA;
  ///done for potential

  Return_t evaluate(ParticleSet& P);
  void addEnergy(MCWalkerConfiguration &W,
                 std::vector<RealType> &LocalEnergy);

  inline Return_t evaluate(ParticleSet& P, std::vector<NonLocalData>& Txy)
  {
    return evaluate(P);
  }

  void addObservables(PropertySetType& plist) { }
  void addObservables(PropertySetType& plist,BufferType& olist);
  void registerCollectables(std::vector<observable_helper*>& h5desc, hid_t gid) const ;
  void setObservables(PropertySetType& plist);
  void setParticlePropertyList(PropertySetType& plist, int offset);
  bool put(xmlNodePtr cur);
  bool get(std::ostream& os) const;
  QMCHamiltonianBase* makeClone(ParticleSet& qp, TrialWaveFunction& psi);

  inline int getGridIndex(int i, int j, int k) const
  {
    return myIndex+k+NumGrids[2]*(j+NumGrids[1]*i);
  }

  inline int getGridIndexPotential(int i, int j, int k) const
  {
    return potentialIndex+k+NumGrids[2]*(j+NumGrids[1]*i);
  }


private:
  ///true if any direction of a supercell is periodic
  bool Periodic;
  ///normalization factor
  RealType Norm;
  ///number of grids
  TinyVector<int,OHMMS_DIM+1> NumGrids;
  ///bin size
  TinyVector<RealType,OHMMS_DIM> Delta;
  ///inverse
  TinyVector<RealType,OHMMS_DIM> DeltaInv;
  ///scaling factor for conversion
  TinyVector<RealType,OHMMS_DIM> ScaleFactor;
  ///lower bound
  TinyVector<RealType, OHMMS_DIM> density_min;
  ///upper bound
  TinyVector<RealType, OHMMS_DIM> density_max;
  ///name of the density data
  std::string prefix;
  ///density
  //Array<RealType,OHMMS_DIM> density, Vavg;
  /** resize the internal data
   *
   * The argument list is not completed
   */
  void resize();
};

}
#endif

