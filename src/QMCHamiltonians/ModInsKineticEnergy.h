#ifndef QMCPLUSPLUS_MODINSKINETICENERGY_H
#define QMCPLUSPLUS_MODINSKINETICENERGY_H
#include "Particle/ParticleSet.h"
#include "Particle/WalkerSetRef.h"
#include "QMCHamiltonians/QMCHamiltonianBase.h"
#include "QMCHamiltonians/MomentumDistribution.h"
#include "QMCWaveFunctions/TrialWaveFunction.h"
#include "OhmmsData/ParameterSet.h"

namespace qmcplusplus {
  
  template<class MomentumDistribution>
  class ModInsKineticEnergy : public QMCHamiltonianBase {
  private:
    // Should results be accumulated?
    bool Accumulate;
    // Total Number of displacements in the current Momentum Distribution
    IndexType NumSamples;
    // Number of displacements used in each call to MomentumDistribution.update()
    IndexType NumDispl;
    // One over the Mass of the particle
    RealType OneOverM;
    // Dispersion Relation holds a list of values giving KE(k)
    Vector<RealType> DispRel;
    // Momentum Distribution is used to calculate n(k) for a given 
    // configuration or group of configurations
    MomentumDistribution MomDist;
    // Reference to the wave function
    TrialWaveFunction& WaveFcn;
    // Particle Set Reference For use in cloning
    ParticleSet& PtclSet;
    
  public:
    /** constructor
     *
     * Kinetic energy operators need to be re-evaluated during optimization.
     */
    ModInsKineticEnergy(ParticleSet& p, TrialWaveFunction& psi, 
                        const Vector<RealType>& DispersRelat,
                        const Vector<IndexType>& DimSizes) :
    Accumulate(false), NumSamples(0), NumDispl(1), OneOverM(1.0), DispRel(DispersRelat), 
    MomDist(DimSizes), WaveFcn(psi), PtclSet(p) {
      UpdateMode.set(OPTIMIZABLE,1);
    }
    
    // copy constructor
    ModInsKineticEnergy(const ModInsKineticEnergy& rhs) : Accumulate(rhs.Accumulate),
    NumSamples(rhs.NumSamples), OneOverM(rhs.OneOverM), DispRel(rhs.DispRel),
    MomDist(rhs.MomDist), WaveFcn(rhs.WaveFcn), PtclSet(rhs.PtclSet) {
      UpdateMode.set(OPTIMIZABLE,1);
    }
    
    ///destructor
    ~ModInsKineticEnergy() { }
    
    inline void SetAccumulate(bool newValue) { Accumulate = newValue; }
    
    void resetTargetParticleSet(ParticleSet& P) { }
 
    inline Return_t evaluate(ParticleSet& P) {
      if(Accumulate) {
        MomDist.updateDistribution(P, WaveFcn, NumDispl);
        Value = 0.0;
      } else {
        MomDist.resetDistribution();
        MomDist.updateDistribution(P, WaveFcn, NumDispl);
        Vector<RealType>& NofK(MomDist.getNofK());
        RealType retVal = 0.0;
        for (IndexType i = 0; i < NofK.size(); i++) {
          retVal += NofK[i] * DispRel[i];
        }
        Value = OneOverM * retVal;
      }
      return Value;
    }
    
    Return_t getEnsembleAverage() {
      Vector<RealType>& NofK(MomDist.getNofK());
      MomDist.resetDistribution();
      RealType retVal = 0.0;
      for (IndexType i = 0; i < NofK.size(); i++) {
        retVal += NofK[i] * DispRel[i];
      }
//      std::cout << retVal * OneOverM << std::endl;
      return OneOverM * retVal;
    }
    
    /** implements the virtual function.
     * 
     * Needs to output all of the constructor arguments as well as template arguments
     */
    bool put(xmlNodePtr cur) {
      ParameterSet newParam;
      newParam.add(NumDispl, "numSamples", "int");
      newParam.add(Accumulate, "accumulate", "bool");
      newParam.add(OneOverM, "invMass", "double");
      newParam.put(cur);
      return true;
    }
    
    bool get(std::ostream& os) const {
      os << "Model Insulator Kinetic Energy";
      return true;
    }
    
    QMCHamiltonianBase* clone() {
      return new ModInsKineticEnergy<MomentumDistribution>(*this);
    }
  };

}

#endif
