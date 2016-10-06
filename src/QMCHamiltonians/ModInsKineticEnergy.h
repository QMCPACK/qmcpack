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
    
    
#ifndef QMCPLUSPLUS_MODINSKINETICENERGY_H
#define QMCPLUSPLUS_MODINSKINETICENERGY_H
#include "Particle/ParticleSet.h"
#include "Particle/WalkerSetRef.h"
#include "QMCHamiltonians/QMCHamiltonianBase.h"
#include "QMCHamiltonians/MomentumDistribution.h"
#include "QMCWaveFunctions/TrialWaveFunction.h"
#include "OhmmsData/ParameterSet.h"

namespace qmcplusplus
{

class ModInsKineticEnergy : public QMCHamiltonianBase
{
private:
  // Should results be accumulated?
  bool Accumulate;
  // Total Number of displacements in the current Momentum Distribution
  IndexType NumSamples;
  // Number of Cycles through the grid used in each updateDistribution call
  IndexType NumCycles;
  // One over the Mass of the particle
  RealType OneOverM;
  // Dispersion Relation holds a list of values giving KE(k)
  Vector<RealType> DispRel;
  // Momentum Distribution is used to calculate n(k) for a given
  // configuration or group of configurations
  MomDistBase* MomDistP;
  // Reference to the wave function
  TrialWaveFunction& WaveFcn;

public:
  /** constructor
   *
   * Kinetic energy operators need to be re-evaluated during optimization.
   */
  ModInsKineticEnergy(TrialWaveFunction& psi, const Vector<RealType>& DispersRelat,
                      const MomDistBase* MomPtr) :
    Accumulate(false), NumSamples(0), NumCycles(1), OneOverM(1.0), DispRel(DispersRelat),
    WaveFcn(psi)
  {
    UpdateMode.set(OPTIMIZABLE,1);
    MomDistP = MomPtr->clone();
  }

  // copy constructor
  ModInsKineticEnergy(const ModInsKineticEnergy& rhs) : Accumulate(rhs.Accumulate),
    NumSamples(rhs.NumSamples), NumCycles(rhs.NumCycles), OneOverM(rhs.OneOverM),
    DispRel(rhs.DispRel), WaveFcn(rhs.WaveFcn)
  {
    UpdateMode.set(OPTIMIZABLE,1);
    MomDistP = rhs.MomDistP->clone();
  }

  ///destructor
  ~ModInsKineticEnergy()
  {
    if (MomDistP)
      delete MomDistP;
  }

  inline void SetAccumulate(bool newValue)
  {
    Accumulate = newValue;
  }

  void resetTargetParticleSet(ParticleSet& P) { }

  inline Return_t evaluate(ParticleSet& P)
  {
    if(Accumulate)
    {
      MomDistP->updateDistribution(P, WaveFcn, NumCycles);
      Value = 0.0;
    }
    else
    {
      MomDistP->resetDistribution();
      MomDistP->updateDistribution(P, WaveFcn, NumCycles);
      Vector<RealType>& NofK(MomDistP->getNofK());
      RealType retVal = 0.0;
      for (IndexType i = 0; i < NofK.size(); i++)
      {
        retVal += NofK[i] * DispRel[i];
      }
      Value = OneOverM * retVal;
    }
    return Value;
  }

  inline Return_t evaluate(ParticleSet& P, std::vector<NonLocalData>& Txy)
  {
    return evaluate(P);
  }

  Return_t getEnsembleAverage()
  {
    Vector<RealType>& NofK(MomDistP->getNofK());
    RealType retVal = 0.0;
    for (IndexType i = 0; i < NofK.size(); i++)
    {
      retVal += NofK[i] * DispRel[i];
//        std::cout << i << "    " << NofK[i] << "     " << DispRel[i] << "    " << retVal << std::endl;
    }
//      std::cout << std::endl;
    MomDistP->resetDistribution();
    return OneOverM * retVal;
  }

  /** implements the virtual function.
   *
   * Needs to output all of the constructor arguments as well as template arguments
   */
  bool put(xmlNodePtr cur)
  {
    ParameterSet newParam;
    newParam.add(NumCycles, "numSamples", "int");
    newParam.add(Accumulate, "accumulate", "bool");
    newParam.add(OneOverM, "invMass", "double");
    newParam.put(cur);
    return true;
  }

  bool get(std::ostream& os) const
  {
    os << "Model Insulator Kinetic Energy";
    return true;
  }

  QMCHamiltonianBase* clone()
  {
    return new ModInsKineticEnergy(*this);
  }
};

}

#endif
