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

namespace qmcplusplus {

  // particle choice policies
  class RandomChoicePolicy : public QMCTraits {
  private:
    const IndexType NumParticles;
  public:
    RandomChoicePolicy(ParticleSet& p) : NumParticles(p.getTotalNum()) { ; }
    RandomChoicePolicy(const RandomChoicePolicy& arg) : NumParticles(arg.NumParticles) { ; }
    inline void NewWalker() { ; }
    inline IndexType operator() () {
      return static_cast<IndexType>(Random() * NumParticles);
    }
  };

  class RandomChoicePerWalkerPolicy : public QMCTraits {
  private:
    bool IsCurrent;
    IndexType CurrentChoice;
    const IndexType NumParticles;
  public:
    RandomChoicePerWalkerPolicy(ParticleSet& p) : IsCurrent(false), NumParticles(p.getTotalNum()) { ; }
    RandomChoicePerWalkerPolicy(const RandomChoicePerWalkerPolicy& rhs) :
      IsCurrent(rhs.IsCurrent), CurrentChoice(rhs.CurrentChoice), NumParticles(rhs.NumParticles) { ; }
    inline void NewWalker() { IsCurrent = false; }
    inline IndexType operator() () {
      if (!IsCurrent) {
        IsCurrent = true;
        CurrentChoice = static_cast<IndexType>(Random() * NumParticles);
      } 
      return CurrentChoice;
    }
  };
  
  class StaticChoicePolicy : public QMCTraits {
  private:
    IndexType CurrentChoice;
  public:
    StaticChoicePolicy(ParticleSet& p) { 
      CurrentChoice = static_cast<IndexType>(Random() * p.getTotalNum()); 
    }
    StaticChoicePolicy(const StaticChoicePolicy& rhs) : CurrentChoice(rhs.CurrentChoice) { ; }
    inline void NewWalker() { ; }
    inline IndexType operator() () {
      return CurrentChoice;
    }
  };
  
  class MomDistBase : public QMCTraits {
  protected:
    IndexType totalNumSamples;
    TinyVector<IndexType, 3> NumPts;
    PosType InvSpacing;
    FFTAbleVectorBase<ComplexType>* NofK;
    Vector<RealType> MomDist;

    template<IndexType DistrDims>
    void placeInBin(const PosType& FracDisp, ValueType RatioOfWvFcn) { 
      TinyVector<IndexType, DistrDims> locInDistribution;
      for (IndexType i = 0; i < DistrDims; i++) {
	locInDistribution[i] = IndexType(std::floor(FracDisp[i] * NumPts[i]));
      }
      IndexType IndexInNofK = locInDistribution[DistrDims-1];
      for (IndexType i = 1; i < DistrDims; i++) {
	IndexInNofK *= NumPts[DistrDims-i];
	IndexInNofK += locInDistribution[DistrDims-i-1];
      }
      (*NofK)(IndexInNofK) += RatioOfWvFcn;
      totalNumSamples++;
    }
     
    template<IndexType DisplDims, IndexType DistrDims> 
    void updateDistHelper(ParticleSet& PtclSet, TrialWaveFunction& Psi, IndexType PtclToDisplace) {
      PosType dr;
      PosType fracDisp;
      for (IndexType i = 0; i < DisplDims; i++) {
	fracDisp[i] = Random();
	dr += fracDisp[i] * PtclSet.Lattice.a(i);
      }
      PtclSet.makeMove(PtclToDisplace, dr);
      placeInBin<DistrDims>(fracDisp, Psi.ratio(PtclSet, PtclToDisplace));
      PtclSet.rejectMove(PtclToDisplace);
    }
  public:
    // It is expected that all of the values will be initialized in the derived classes!
    MomDistBase() { ; }
    ~MomDistBase() {
      delete NofK;
    }
    MomDistBase(const MomDistBase& arg) : totalNumSamples(arg.totalNumSamples),
    InvSpacing(arg.InvSpacing), MomDist(arg.MomDist) { 
       NofK = arg.NofK->clone();
    }
    
    inline void resetDistribution() { *NofK *= 0.0; totalNumSamples = 0; }
     
    virtual void updateDistribution(ParticleSet& PtclSet, TrialWaveFunction& Psi, IndexType NumSamples) = 0;
     
    Vector<RealType>& getNofK() {
      NofK->setForwardNorm(1.0 / (NofK->size() * totalNumSamples));
      NofK->transformForward();
      for (IndexType i = 0; i < NofK->size(); i++) {
        MomDist[i] = (*NofK)[i].real();
      }
      return MomDist;
    }
  };	  
   
  // really ugly, but works for now and shows the code that needs to
  // be put some other place eventually
  template<> void MomDistBase::updateDistHelper<3,1>(ParticleSet& PtclSet, TrialWaveFunction& Psi, IndexType PtclToDisplace) {
    PosType dr;
    for (IndexType i = 0; i < 3; i++) {
      dr += Random() * PtclSet.Lattice.a(i);
    }
    PtclSet.makeMove(PtclToDisplace, dr);
    RealType DispVecLength = std::sqrt(dot(dr, dr));
    IndexType IndexInNofK = int(std::floor(NumPts[0] * InvSpacing[0] * DispVecLength));
    if (IndexInNofK < NumPts[0]) {
      (*NofK)[IndexInNofK] += Psi.ratio(PtclSet, PtclToDisplace);
       totalNumSamples++;
    }
    PtclSet.rejectMove(PtclToDisplace);
  }

  template<class PtclChoicePolicy>
  class ThreeDimMomDist : public MomDistBase {
  private:
    PtclChoicePolicy pcp;
  public:
    ThreeDimMomDist(ParticleSet& PtclSet, const Vector<IndexType>& DimSizes) : pcp(PtclSet) {
      totalNumSamples = 0;
      for (IndexType i = 0; i < 3; i++) {
        NumPts[i] = DimSizes[i];
        InvSpacing = 1.0 / std::sqrt(dot(PtclSet.Lattice.a(i), PtclSet.Lattice.a(i)));
      }
      NofK = new FFTAbleVector<3, ComplexType, FFTWEngine>(NumPts);
      MomDist.resize(NofK->size());
    }
    ThreeDimMomDist(const ThreeDimMomDist& arg) : MomDistBase(arg), pcp(arg.pcp) { ; }
    void updateDistribution(ParticleSet& PtclSet, TrialWaveFunction& Psi, IndexType NumSamples) {
      IndexType targetNumSamples = totalNumSamples + NumSamples;
      while(totalNumSamples < targetNumSamples) {
        IndexType TargetParticle = pcp();
        updateDistHelper<3,3>(PtclSet, Psi, TargetParticle);
      }
    }
  };

  template<class PtclChoicePolicy>
  class OneDimMomDist : public MomDistBase {
  private:
    PtclChoicePolicy pcp;
  public:
    OneDimMomDist(ParticleSet& PtclSet, const Vector<IndexType>& DimSizes) : pcp(PtclSet) {
      totalNumSamples = 0;
      for (IndexType i = 0; i < 1; i++) {
        NumPts[i] = DimSizes[i];
        InvSpacing = 1.0 / std::sqrt(dot(PtclSet.Lattice.a(i), PtclSet.Lattice.a(i)));
      }
      NofK = new FFTAbleVector<1, ComplexType, FFTWEngine>(NumPts);
      MomDist.resize(NofK->size());
    }
    OneDimMomDist(const OneDimMomDist& arg) : MomDistBase(arg), pcp(arg.pcp) { ; }
    void updateDistribution(ParticleSet& PtclSet, TrialWaveFunction& Psi, IndexType NumSamples) {
      IndexType targetNumSamples = totalNumSamples + NumSamples;
      while(totalNumSamples < targetNumSamples) {
        IndexType TargetParticle = pcp();
        updateDistHelper<1,1>(PtclSet, Psi, TargetParticle);
      }
    }
  };

  template<class PtclChoicePolicy>
  class AveragedOneDimMomDist : public MomDistBase {
  private:
    PtclChoicePolicy pcp;
  public:
    AveragedOneDimMomDist(ParticleSet& PtclSet, const Vector<IndexType>& DimSizes) : pcp(PtclSet) {
      totalNumSamples = 0;
      for (IndexType i = 0; i < 1; i++) {
        NumPts[i] = DimSizes[i];
        InvSpacing = 1.0 / std::sqrt(dot(PtclSet.Lattice.a(i), PtclSet.Lattice.a(i)));
      }
      NofK = new FFTAbleVector<1, ComplexType, FFTWEngine>(NumPts);
      MomDist.resize(NofK->size());
    }
    AveragedOneDimMomDist(const AveragedOneDimMomDist& arg) : MomDistBase(arg), pcp(arg.pcp) { ; }
    void updateDistribution(ParticleSet& PtclSet, TrialWaveFunction& Psi, IndexType NumSamples) {
      IndexType targetNumSamples = totalNumSamples + NumSamples;
      while(totalNumSamples < targetNumSamples) {
        IndexType TargetParticle = pcp();
        updateDistHelper<3,1>(PtclSet, Psi, TargetParticle);
      }
    }
  };

}
           
#endif
