#ifndef QMCPLUSPLUS_COUNTINGJASTROWORBITAL_H
#define QMCPLUSPLUS_COUNTINGJASTROWORBITAL_H
#include "QMCWaveFunctions/OrbitalBase.h"
#include "QMCWaveFunctions/Jastrow/CountingRegion.h"
#include "QMCWaveFunctions/Jastrow/CountingFunctors.h"

namespace qmcplusplus
{
class CountingJastrowOrbital : public OrbitalBase
{
protected:
  // number of electrons
  int num_els;
  // number of counting regions
  int num_regions;
  // debug print flag
  bool debug;
  int debug_seqlen;
  int debug_period;

  // Jastrow linear coefficients
  std::vector<RealType> F;
  std::vector<RealType> G;
  // Counting Regions
  CountingRegionBase* C;

  // Optimization Flags
  bool opt_F;
  bool opt_G;
  bool opt_C;

  // flag for using normalized counting regions
  bool C_norm;

  // Jastrow intermediate Matrix-vector products 
  std::vector<RealType> FCsum;
  std::vector<GradType> FCgrad;
  std::vector<RealType> FClap;
  
  // grad dot grad and laplacian sums for evaluateDerivatives
  std::vector<RealType> FCggsum;
  std::vector<RealType> FClapsum;

  // Jastrow intermediate Matrix-vector products at proposed position
  std::vector<RealType> FCsum_t;
  std::vector<GradType> FCgrad_t;
  std::vector<RealType> FClap_t;
  
  // Jastrow exponent values and gradients (by particle index)
  RealType Jval; 
  std::vector<GradType> Jgrad;
  std::vector<RealType> Jlap;

  // Jastrow exponent values and gradients at proposed position
  RealType Jval_t;
  std::vector<GradType> Jgrad_t;
  std::vector<RealType> Jlap_t;

  // containers for counting function derivative quantities
  std::vector<RealType> dCsum;
  std::vector<RealType> dCsum_t;
  std::vector<RealType> dCggsum;
  std::vector<RealType> dClapsum;
  std::vector<RealType> dCFCggsum;
  std::vector<int> dCindex;


  
  // first array index for opt_index, opt_id 
  enum opt_var { OPT_F, OPT_G, NUM_OPT_VAR };
  // vectors to store indices and names of active optimizable parameters  
  std::array<std::vector<int>,NUM_OPT_VAR> opt_index; 
  std::array<std::vector<std::string>,NUM_OPT_VAR> opt_id;

public:
  CountingJastrowOrbital(ParticleSet& els);
  ~CountingJastrowOrbital();
  bool put(xmlNodePtr cur);
  void initialize();
  void reportStatus(std::ostream& os);

  // internally used functions
  void evaluateExponents(ParticleSet& P);
  void evaluateTempExponents(ParticleSet& P, int iat);

  // print helper functions
  void evaluateExponents_print(std::ostream& os, ParticleSet& P);
  void evaluateTempExponents_print(std::ostream& os, ParticleSet& P, int iat);


  // === VMC calls ===  
  // called at the beginning of every VMC run
  void resetTargetParticleSet(ParticleSet& P);
  inline RealType registerData(ParticleSet& P, PooledData<RealType>& buf);
  // called in VMC runs before each step
  inline void copyFromBuffer(ParticleSet& P,PooledData<RealType>& buf);
  // called in VMC runs in advanceWalkers: return Gradient of particle at iat
  GradType evalGrad(ParticleSet& P, int iat);
  // called in VMC runs in advanceWalkers: calculate wfn and ratio at particle position iat
  ValueType ratioGrad(ParticleSet& P, int iat, GradType& grad_iat);
  // called in VMC runs after each substep if step is accepted
  void acceptMove(ParticleSet& P, int iat);
  // called in VMC runs after each substep if step is rejected
  void restore(int iat);
  // called in VMC runs after each step  
  inline RealType updateBuffer(ParticleSet& P, PooledData<RealType>& buf, bool fromscratch=false);
  // update not called?
  void update(ParticleSet& P, ParticleSet::ParticleGradient_t& dG, ParticleSet::ParticleLaplacian_t& dL, int iat);
  // called every nBlocksBetweenRecompute: default 0
  void recompute(ParticleSet& P);
  // called after recompute
  inline RealType evaluateLog(ParticleSet& P, PooledData<RealType>& buf);

  // === linear method calls ===
  // called from QMCCostFunctionBase.put: once at the beginning of each qmcDriver
  void checkInVariables(opt_variables_type& active);
  void checkOutVariables(const opt_variables_type& active);
  // called from LMYEngineCost, four times (initial position, three shifts)
  void resetParameters(const opt_variables_type& active);
  // called in evaluateDeltaLog. Primarly linear method check configurations
  RealType evaluateLog(ParticleSet& P, 
                       ParticleSet::ParticleGradient_t& G,
                       ParticleSet::ParticleLaplacian_t& L);
  // called in evaluateValueAndDerivatives
  ValueType ratio(ParticleSet& P, int iat); 
  
  // derivatives of single-particle update for ecp nonlocal quadrature
  void evaluateTempDerivatives(ParticleSet& P, 
                                const opt_variables_type& active, 
                                RealType& psiratio,
                                std::vector<RealType>& dlogpsi_t, 
                                int iel,
                                PosType dr);

  void evaluateDerivatives(ParticleSet& P, 
                           const opt_variables_type& active, 
                           std::vector<RealType>& dlogpsi, 
                           std::vector<RealType>& dhpsioverpsi);

  // === DMC method calls === 
  // called in DMC runs in advanceWalkers
  ValueType ratio(ParticleSet& P, int iat, ParticleSet::ParticleGradient_t& dG,ParticleSet::ParticleLaplacian_t& dL);

  // storageType = 0 ; store everything. Inherits no 
  // void registerDataForDerivatives(ParticleSet& P, PooledData<RealType>& buf, int storageType)

  // Psi->evaluate(P) // not used? required for basis particle updaters
  ValueType evaluate(ParticleSet& P, ParticleSet::ParticleGradient_t& G, ParticleSet::ParticleLaplacian_t& L);

  OrbitalBasePtr makeClone(ParticleSet &P) const;


};

}

#endif
