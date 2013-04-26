#ifndef QMCPLUSPLUS_RPA_JASTROW_H
#define QMCPLUSPLUS_RPA_JASTROW_H

#include "QMCWaveFunctions/OrbitalBase.h"
#include "LongRange/LRHandlerBase.h"
#include "QMCWaveFunctions/Jastrow/SplineFunctors.h"
#include "QMCWaveFunctions/Jastrow/TwoBodyJastrowOrbital.h"
#include "QMCWaveFunctions/Jastrow/LRBreakupUtilities.h"
#include "QMCWaveFunctions/Jastrow/LRTwoBodyJastrow.h"


namespace qmcplusplus
{

/** JastrowBuilder using RPA functor
 *  Modification of RPAJastrow
 *
 */
struct RPAJastrow: public OrbitalBase
{
  typedef LRHandlerBase HandlerType;
  typedef CubicBspline<RealType,LINEAR_1DGRID,FIRSTDERIV_CONSTRAINTS> SplineEngineType;
  typedef CubicSplineSingle<RealType,SplineEngineType> FuncType;
  typedef LinearGrid<RealType> GridType;

  RPAJastrow(ParticleSet& target, bool is_manager);

  ~RPAJastrow();

  bool put(xmlNodePtr cur);

  void buildOrbital(const string& name, const string& UL
                    , const string& US, const string& RF, RealType R, RealType K);

  void makeShortRange();
  void makeLongRange();

  void setHandler(HandlerType* Handler)
  {
    myHandler=Handler;
  };

  /** check out optimizable variables
    */
  void checkOutVariables(const opt_variables_type& o);

  /** check in an optimizable parameter
        * @param o a super set of optimizable variables
    */
  void checkInVariables(opt_variables_type& o);

  /** print the state, e.g., optimizables */
  void reportStatus(ostream& os);

  /** reset the parameters during optimizations
    */
  void resetParameters(const opt_variables_type& active);

  void resetTargetParticleSet(ParticleSet& P);

  ValueType evaluate(ParticleSet& P,
                     ParticleSet::ParticleGradient_t& G,
                     ParticleSet::ParticleLaplacian_t& L)
  {
    return std::exp(evaluateLog(P,G,L));
  }

  RealType evaluateLog(ParticleSet& P,
                       ParticleSet::ParticleGradient_t& G, ParticleSet::ParticleLaplacian_t& L);

  ValueType ratio(ParticleSet& P, int iat,
                  ParticleSet::ParticleGradient_t& dG,
                  ParticleSet::ParticleLaplacian_t& dL);

  ValueType ratio(ParticleSet& P, int iat);

  void acceptMove(ParticleSet& P, int iat);

  void restore(int iat);

  void update(ParticleSet& P,
              ParticleSet::ParticleGradient_t& dG,
              ParticleSet::ParticleLaplacian_t& dL,
              int iat);

  RealType registerData(ParticleSet& P, BufferType& buf);

  RealType updateBuffer(ParticleSet& P, BufferType& buf, bool fromscratch=false);

  void copyFromBuffer(ParticleSet& P, BufferType& buf);

  RealType evaluateLog(ParticleSet& P,BufferType& buf);

  OrbitalBase* makeClone(ParticleSet& tqp) const;

private:

  bool IsManager;
  bool IgnoreSpin;
  bool DropLongRange;
  bool DropShortRange;
  RealType Rs;
  RealType Kc;
  RealType Rcut;
  string ID_Rs;
  string rpafunc;
  string MyName;

  ///@}
  /** main handler
   */
  HandlerType* myHandler;
  ///object to handle the long-range part
  LRTwoBodyJastrow* LongRangeRPA;
  ///@{objects to handle the short-range part
  ///two-body Jastrow function
  OrbitalBase* ShortRangeRPA;
  ///numerical function owned by ShortRangeRPA
  FuncType* nfunc;
  GridType* myGrid;
  ///adaptor function to initialize nfunc
  ShortRangePartAdapter<RealType>* SRA;
  ParticleSet& targetPtcl;
  ///A list of OrbitalBase*
  vector<OrbitalBase*> Psi;
};
}
#endif
