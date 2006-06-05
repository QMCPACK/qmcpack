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
  public:
    void NewWalker() { ; }
    IndexType operator() (const ParticleSet& p) {
      return static_cast<IndexType>(Random() * p.getTotalNum());
    }
  };

  class RandomChoicePerWalkerPolicy : public QMCTraits {
  private:
    bool IsCurrent;
    IndexType CurrentChoice;
  public:
    RandomChoicePerWalkerPolicy() : IsCurrent(false) { ; }
    RandomChoicePerWalkerPolicy(const RandomChoicePerWalkerPolicy& rhs) :
      IsCurrent(rhs.IsCurrent), CurrentChoice(rhs.CurrentChoice) { ; }
    void NewWalker() { IsCurrent = false; }
    IndexType operator() (const ParticleSet& p) {
      if (!IsCurrent) {
        IsCurrent = true;
        CurrentChoice = static_cast<IndexType>(Random() * p.getTotalNum());
      } 
      return CurrentChoice;
    }
  };
  
  class StaticChoicePolicy : public QMCTraits {
  private:
    bool IsCurrent;
    IndexType CurrentChoice;
  public:
    StaticChoicePolicy() : IsCurrent(false) { ; }
    StaticChoicePolicy(const StaticChoicePolicy& rhs) :
      IsCurrent(rhs.IsCurrent), CurrentChoice(rhs.CurrentChoice) { ; }
    void NewWalker() { ; }
    IndexType operator() (const ParticleSet& p) {
      if (!IsCurrent) {
        IsCurrent = true;
        CurrentChoice = static_cast<IndexType>(Random() * p.getTotalNum());
      } 
      return CurrentChoice;
    }
  };
  
  // Momentum Distributions
  template<class PtclChoicePolicy>
  class ThreeDimMomDist : public QMCTraits {
  private:
    IndexType totalNumSamples;
    IndexType nx;
    IndexType ny; 
    IndexType nz;
    PtclChoicePolicy pcp;
    FFTAbleVector<3, ComplexType, FFTWEngine>* NofK;
    void placeInBin(const ParticleSet& PtclSet, const PosType& FracDispVector, ValueType RatioOfWvFcn) {
      RealType xSpacing = std::sqrt(dot(PtclSet.Lattice.a(0), PtclSet.Lattice.a(0))) / nx;
      RealType ySpacing = std::sqrt(dot(PtclSet.Lattice.a(1), PtclSet.Lattice.a(1))) / ny;
      RealType zSpacing = std::sqrt(dot(PtclSet.Lattice.a(2), PtclSet.Lattice.a(2))) / nz;
      IndexType locx = int(std::floor(FracDispVector[0] / xSpacing));
      IndexType locy = int(std::floor(FracDispVector[1] / ySpacing));
      IndexType locz = int(std::floor(FracDispVector[2] / zSpacing));
      (*NofK)(locz + locy * nz + locx * ny * nz) += RatioOfWvFcn;
    }
  public:
    ThreeDimMomDist(const Vector<IndexType>& DimSizes) : totalNumSamples(0), nx(0), ny(0), nz(0) {
      nx = DimSizes[0] / 2;
      ny = DimSizes[1] / 2;
      nz = DimSizes[2] / 2;
      NofK = new FFTAbleVector<3, ComplexType, FFTWEngine>(nx, ny, nz);
    }
    
    ThreeDimMomDist(const ThreeDimMomDist& rhs) : totalNumSamples(rhs.totalNumSamples),
    nx(rhs.nx), ny(rhs.ny), nz(rhs.nz), pcp(rhs.pcp) {
      NofK = new FFTAbleVector<3, ComplexType, FFTWEngine>(*(rhs.NofK));
    }
    
    ~ThreeDimMomDist() {
      delete NofK;
    }
    
    inline void resetDistribution() { *NofK *= 0.0; totalNumSamples = 0; }
    
    void updateDistribution(ParticleSet& PtclSet, TrialWaveFunction& Psi, IndexType NumSamples) {
      PosType dr;
      PosType fracDisp;
      // Tell the particle choice policy that this is for a new walker
      pcp.NewWalker();
      // for NumSamplesTimes
      totalNumSamples += NumSamples;
      for (IndexType sample = 0; sample < NumSamples; sample++) {
        // choose random displacement in the box
        dr[0] = dr[1] = dr[2] = 0.0;
        for (IndexType i = 0; i < 3; i++) {
          fracDisp[i] = Random();
          dr += fracDisp[i] * PtclSet.Lattice.a(i);
        }
        IndexType particleIndex = pcp(PtclSet);
        PtclSet.makeMove(particleIndex, dr);
        placeInBin(PtclSet, fracDisp, Psi.ratio(PtclSet, particleIndex));
        PtclSet.rejectMove(particleIndex);
      }
    }
    
    Vector<ValueType> getNofK() {
      NofK->setForwardNorm(1.0 / NofK->size());
      NofK->transformForward();
      Vector<ValueType> temp;
      temp.resize(NofK->size());
      for (IndexType i = 0; i < NofK->size(); i++) {
        temp[i] = (*NofK)[i].real();
      }
      return temp;
    }
  };

  template<class PtclChoicePolicy>
  class OneDimMomDist : public QMCTraits {
  private:
    IndexType totalNumSamples;
    IndexType nx;
    PtclChoicePolicy pcp;
    FFTAbleVector<1, ComplexType, FFTWEngine>* NofK;
    void placeInBin(const ParticleSet& PtclSet, const PosType& DispVector, ValueType RatioOfWvFcn) {
      RealType DispVecLength = std::sqrt(dot(DispVector, DispVector));
      RealType Spacing = std::sqrt(dot(PtclSet.Lattice.a(0), PtclSet.Lattice.a(0))) / nx;
      IndexType xloc = int(std::floor(DispVecLength / Spacing));
      (*NofK)[xloc] += RatioOfWvFcn;
    } 
  public:
    OneDimMomDist(const Vector<IndexType>& DimSizes) : totalNumSamples(0), nx(0) {
      nx = DimSizes[0] / 2;
      NofK = new FFTAbleVector<1, ComplexType, FFTWEngine>(nx);
    }

    OneDimMomDist(const OneDimMomDist& rhs) : totalNumSamples(rhs.totalNumSamples), 
    nx(rhs.nx), pcp(rhs.pcp) {
      NofK = new FFTAbleVector<1, ComplexType, FFTWEngine>(*(rhs.NofK));
    }
    
    ~OneDimMomDist() {
      delete NofK;
    }
    
    inline void resetDistribution() { *NofK *= 0.0; totalNumSamples = 0; }
    
    void updateDistribution(ParticleSet& PtclSet, TrialWaveFunction& Psi, IndexType NumSamples) {
      PosType dr;
      // Tell the particle choice policy that this is for a new walker
      pcp.NewWalker();
      // for NumSamplesTimes
      totalNumSamples += NumSamples;
      for (IndexType sample = 0; sample < NumSamples; sample++) {
        // choose random displacement along line
        dr = Random() * PtclSet.Lattice.a(0);
        IndexType particleIndex = pcp(PtclSet);
        PtclSet.makeMove(particleIndex, dr);
        placeInBin(PtclSet, dr, Psi.ratio(PtclSet, particleIndex));
        PtclSet.rejectMove(particleIndex);
      }
    }
    
    Vector<ValueType> getNofK() {
      NofK->setForwardNorm(1.0 / NofK->size());
      NofK->transformForward();
      Vector<ValueType> temp;
      temp.resize(NofK->size());
      for (IndexType i = 0; i < NofK->size(); i++) {
        temp[i] = (*NofK)[i].real();
      }
      return temp;
    }
  };
  
  
  template<class PtclChoicePolicy>
  class AveragedOneDimMomDist : public QMCTraits {
  private:
    IndexType totalNumSamples;
    IndexType nx;
    PtclChoicePolicy pcp;
    FFTAbleVector<1, ComplexType, FFTWEngine>* NofK;
    void placeInBin(const ParticleSet& PtclSet, const PosType& DispVector, ValueType RatioOfWvFcn) {
      RealType DispVecLength = std::sqrt(dot(DispVector, DispVector));
      RealType Spacing = std::sqrt(dot(PtclSet.Lattice.a(0), PtclSet.Lattice.a(0))) / nx;
      IndexType xloc = int(std::floor(DispVecLength / Spacing));
      if (xloc < nx) (*NofK)[xloc] += RatioOfWvFcn;
    } 
  public:
    AveragedOneDimMomDist(const Vector<IndexType>& DimSizes) : totalNumSamples(0), nx(0) {
      nx = DimSizes[0] / 2;
      NofK = new FFTAbleVector<1, ComplexType, FFTWEngine>(nx);
    }

    AveragedOneDimMomDist(const AveragedOneDimMomDist& rhs) : totalNumSamples(rhs.totalNumSamples), 
    nx(rhs.nx), pcp(rhs.pcp) {
      NofK = new FFTAbleVector<1, ComplexType, FFTWEngine>(*(rhs.NofK));
    }
    
    ~AveragedOneDimMomDist() {
      delete NofK;
    }
    
    inline void resetDistribution() { *NofK *= 0.0; totalNumSamples = 0; }
    
    void updateDistribution(ParticleSet& PtclSet, TrialWaveFunction& Psi, IndexType NumSamples) {
      PosType dr;
      // Tell the particle choice policy that this is for a new walker
      pcp.NewWalker();
      // for NumSamplesTimes
      totalNumSamples += NumSamples;
      for (IndexType sample = 0; sample < NumSamples; sample++) {
        // choose random displacement 
        dr[0] = dr[1] = dr[2] = 0;
        for (IndexType i = 0; i < 3; i++) {
          dr += Random() * PtclSet.Lattice.a(i);
        }
        IndexType particleIndex = pcp(PtclSet);
        PtclSet.makeMove(particleIndex, dr);
        placeInBin(PtclSet, dr, Psi.ratio(PtclSet, particleIndex));
        PtclSet.rejectMove(particleIndex);
      }
    }
    
    Vector<ValueType> getNofK() {
      NofK->setForwardNorm(1.0 / NofK->size());
      NofK->transformForward();
      Vector<ValueType> temp;
      temp.resize(NofK->size());
      for (IndexType i = 0; i < NofK->size(); i++) {
        temp[i] = (*NofK)[i].real();
      }
      return temp;
    }
  };
  
}
           
#endif
