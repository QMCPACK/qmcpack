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
    
    
#ifndef QMCPLUSPLUS_MOMENTUMDISTRIBUTION_H
#define QMCPLUSPLUS_MOMENTUMDISTRIBUTION_H

#include "Configuration.h"
#include "Particle/ParticleSet.h"
#include "QMCWaveFunctions/TrialWaveFunction.h"
#include "OhmmsPETE/OhmmsVector.h"
#include "OhmmsPETE/TinyVector.h"
#include "Utilities/RandomGenerator.h"
#include "Numerics/FFTAbleVector.h"
#include "Numerics/FFTEngines.h"

namespace qmcplusplus
{

class PtclChoiceBase : public QMCTraits
{
public:
  PtclChoiceBase()
  {
    ;
  }
  PtclChoiceBase(const PtclChoiceBase& arg)
  {
    ;
  }
  virtual void NewWalker() = 0;
  virtual IndexType operator() () = 0;
  virtual PtclChoiceBase* clone() const = 0;
  virtual ~PtclChoiceBase()
  {
    ;
  }
};

// particle choice policies
class RandomChoice : public PtclChoiceBase
{
private:
  const IndexType NumParticles;
public:
  RandomChoice(ParticleSet& p) : NumParticles(p.getTotalNum())
  {
    ;
  }
  RandomChoice(const RandomChoice& arg) : NumParticles(arg.NumParticles)
  {
    ;
  }
  ~RandomChoice()
  {
    ;
  }
  inline void NewWalker()
  {
    ;
  }
  inline IndexType operator() ()
  {
    return static_cast<IndexType>(Random() * NumParticles);
  }
  inline RandomChoice* clone() const
  {
    return new RandomChoice(*this);
  };
};

class RandomChoicePerWalker : public PtclChoiceBase
{
private:
  IndexType CurrentChoice;
  const IndexType NumParticles;
public:
  RandomChoicePerWalker(ParticleSet& p) : NumParticles(p.getTotalNum())
  {
    CurrentChoice = static_cast<IndexType>(Random() * NumParticles);
  }
  RandomChoicePerWalker(const RandomChoicePerWalker& rhs) :
    CurrentChoice(rhs.CurrentChoice), NumParticles(rhs.NumParticles)
  {
    ;
  }
  ~RandomChoicePerWalker()
  {
    ;
  }
  inline void NewWalker()
  {
    CurrentChoice = static_cast<IndexType>(Random() * NumParticles);
  }
  inline IndexType operator() ()
  {
    return CurrentChoice;
  }
  inline RandomChoicePerWalker* clone() const
  {
    return new RandomChoicePerWalker(*this);
  };
};

class StaticChoice : public PtclChoiceBase
{
private:
  IndexType CurrentChoice;
public:
  StaticChoice(ParticleSet& p)
  {
    CurrentChoice = static_cast<IndexType>(Random() * p.getTotalNum());
  }
  StaticChoice(const StaticChoice& rhs) : CurrentChoice(rhs.CurrentChoice)
  {
    ;
  }
  ~StaticChoice()
  {
    ;
  }
  inline void NewWalker()
  {
    ;
  }
  inline IndexType operator() ()
  {
    return CurrentChoice;
  }
  inline StaticChoice* clone() const
  {
    return new StaticChoice(*this);
  };
};


class MomDistBase : public QMCTraits
{
protected:
  IndexType NumParticles;
  IndexType totalNumSamples;
  Vector<RealType> MomDist;
  PtclChoiceBase* pcp;

public:
  // It is expected that all of the values will be initialized in the derived classes!
  MomDistBase(PtclChoiceBase* pcb, IndexType NumPart) : NumParticles(NumPart)
  {
    totalNumSamples = 0;
    pcp = pcb->clone();
  }
  virtual ~MomDistBase()
  {
    if (pcp)
      delete pcp;
  }
  MomDistBase(const MomDistBase& arg) : NumParticles(arg.NumParticles), totalNumSamples(arg.totalNumSamples),
    MomDist(arg.MomDist)
  {
    pcp = arg.pcp->clone();
  }

  virtual MomDistBase* clone() const = 0;
  virtual void resetDistribution()
  {
    totalNumSamples = 0;
  }
  virtual void updateDistribution(ParticleSet& PtclSet, TrialWaveFunction& Psi, IndexType NumCycles) = 0;
  virtual Vector<RealType>& getNofK() = 0;
};

template<unsigned dims>
class FFTMomentumDist : public MomDistBase
{
protected:
  TinyVector<IndexType, dims> NumPts;
  FFTAbleVectorBase<ComplexType>* NofK;
  TinyVector<RealType, dims> InvRsLatSize;
  TinyVector<RealType, dims> Spacing;

  // go from a set of integers with respect to NumPts and find the index in
  // the linear structure that holds them
  IndexType getIndex(const TinyVector<IndexType, dims>& LocInDistribution)
  {
    IndexType IndexInNofK = LocInDistribution[0];
    for (IndexType i = 1; i < dims; i++)
    {
      IndexInNofK *= NumPts[i];
      IndexInNofK += LocInDistribution[i];
    }
    return IndexInNofK;
  }

  // These are used when a simple interpolation is needed to put a value that does not
  // fall exactly on the grid onto the grid
  void PlaceCloudInBin(const TinyVector<RealType, dims>& FracDisp, ValueType RatioOfWvFcn);

  // This handles putting a value onto the grid
  inline void placeIntsInBin(const TinyVector<IndexType, dims>& Indexes,
                             ValueType RatioOfWvFcn)
  {
    IndexType index = getIndex(Indexes);
    (*NofK)(index) += RatioOfWvFcn;
  }

  inline void placeInBin(const TinyVector<RealType, dims>& FracDisp, ValueType RatioOfWvFcn)
  {
    TinyVector<IndexType, dims> Indexes;
    for (IndexType i = 0; i < dims; i++)
    {
      Indexes[i] =IndexType(std::floor(FracDisp[i] * NumPts[i]));
    }
    placeIntsInBin(Indexes, RatioOfWvFcn);
  }

  void initialize(const TinyVector<IndexType, dims>& inDimSizes, const TinyVector<RealType, dims>& inInvRsLatSize)
  {
    NumPts = inDimSizes;
    InvRsLatSize = inInvRsLatSize;
    for (IndexType i = 0; i < dims; i++)
    {
      Spacing[i] = 1.0 / (InvRsLatSize[i] * static_cast<RealType>(NumPts[i]));
    }
    NofK = new FFTAbleVector<dims, ComplexType, FFTWEngine>(NumPts);
    MomDist.resize(NofK->size());
  }
  FFTMomentumDist(PtclChoiceBase* pcb, IndexType NumPart) : MomDistBase(pcb, NumPart)
  {
    ;
  }
public:
  ~FFTMomentumDist()
  {
    if (NofK)
      delete NofK;
  }
  FFTMomentumDist(const FFTMomentumDist& arg) : MomDistBase(arg), NumPts(arg.NumPts),
    InvRsLatSize(arg.InvRsLatSize), Spacing(arg.Spacing)
  {
    NofK = arg.NofK->clone();
  }

  void resetDistribution()
  {
    (*NofK) *= 0.0;
    totalNumSamples = 0;
  }

  Vector<RealType>& getNofK()
  {
    NofK->setForwardNorm(static_cast<RealType>(NumParticles) / static_cast<RealType>(totalNumSamples));
    /*
     for(IndexType i = 0; i < NofK->size() / 2; i++) {
       std::cout << (*NofK)[i].real() * static_cast<RealType>(NofK->size()) / static_cast<RealType>(totalNumSamples) << "    ";
     }
     std::cout << std::endl;
    */
    NofK->transformForward();
    for (IndexType i = 0; i < NofK->size(); i++)
    {
      MomDist[i] = (*NofK)[i].real();
      std::cout << MomDist[i] << "    ";
    }
    std::cout << std::endl;
    return MomDist;
  }
};



class ThreeDimMomDist : public FFTMomentumDist<3>
{
public:
  ThreeDimMomDist(ParticleSet& PtclSet, const Vector<IndexType>& DimSizes, PtclChoiceBase* pcb);
  ThreeDimMomDist(const ThreeDimMomDist& arg) : FFTMomentumDist<3>(arg)
  {
    ;
  }
  ~ThreeDimMomDist()
  {
    ;
  }

  void updateDistribution(ParticleSet& PtclSet, TrialWaveFunction& Psi, IndexType NumCycles);
  ThreeDimMomDist* clone() const
  {
    return new ThreeDimMomDist(*this);
  }
};

class OneDimMomDist : public FFTMomentumDist<1>
{
public:
  OneDimMomDist(ParticleSet& PtclSet, const Vector<IndexType>& DimSizes, PtclChoiceBase* pcb);
  OneDimMomDist(const OneDimMomDist& arg) : FFTMomentumDist<1>(arg)
  {
    ;
  }
  ~OneDimMomDist()
  {
    ;
  }

  void updateDistribution(ParticleSet& PtclSet, TrialWaveFunction& Psi, IndexType NumCycles);
  OneDimMomDist* clone() const
  {
    return new OneDimMomDist(*this);
  }
};

class AveragedOneDimMomDist : public FFTMomentumDist<1>
{
public:
  AveragedOneDimMomDist(ParticleSet& PtclSet, const Vector<IndexType>& DimSizes, PtclChoiceBase* pcb);
  AveragedOneDimMomDist(const AveragedOneDimMomDist& arg) : FFTMomentumDist<1>(arg)
  {
    ;
  }
  ~AveragedOneDimMomDist()
  {
    ;
  }

  void updateDistribution(ParticleSet& PtclSet, TrialWaveFunction& Psi, IndexType NumSamples);
  AveragedOneDimMomDist* clone() const
  {
    return new AveragedOneDimMomDist(*this);
  }
};

class DirectMomDist : public MomDistBase
{
protected:
  Vector<PosType> GVectors;
  Vector<ComplexType> NofK;
  DirectMomDist(PtclChoiceBase* pcb, IndexType NumPart) : MomDistBase(pcb, NumPart)
  {
    ;
  }
public:
  DirectMomDist(const DirectMomDist& arg) : MomDistBase(arg), GVectors(arg.GVectors), NofK(arg.NofK)
  {
    ;
  }
  ~DirectMomDist()
  {
    ;
  }
  void resetDistribution()
  {
    totalNumSamples = 0;
    NofK *= 0.0;
  }
  Vector<RealType>& getNofK();
};


class RandomMomDist : public DirectMomDist
{
public:
  RandomMomDist(ParticleSet& PtclSet, const Vector<PosType>& GVectors, PtclChoiceBase* pcb);
  RandomMomDist(const RandomMomDist& arg) : DirectMomDist(arg)
  {
    ;
  }
  ~RandomMomDist()
  {
    ;
  }

  void updateDistribution(ParticleSet& PtclSet, TrialWaveFunction& Psi, IndexType NumSamples);
  RandomMomDist* clone() const
  {
    return new RandomMomDist(*this);
  }
};

}

#endif
