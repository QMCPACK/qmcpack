//////////////////////////////////////////////////////////////////
// (c) Copyright 2003-  by Jeongnim Kim
//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////
//   National Center for Supercomputing Applications &
//   Materials Computation Center
//   University of Illinois, Urbana-Champaign
//   Urbana, IL 61801
//   e-mail: jnkim@ncsa.uiuc.edu
//
// Supported by
//   National Center for Supercomputing Applications, UIUC
//   Materials Computation Center, UIUC
//////////////////////////////////////////////////////////////////
// -*- C++ -*-
/** @file LRTwoBodyJastrow.h
 * @brief Declaration of Long-range TwoBody Jastrow
 *
 * The initial coefficients are evaluated according to RPA.
 * A set of k-dependent coefficients are generated.
 * Optimization should use rpa_k0, ...rpa_kn as the names of the
 * optimizable parameters.
 */
#ifndef QMCPLUSPLUS_LR_RPAJASTROW_H
#define QMCPLUSPLUS_LR_RPAJASTROW_H

#include "QMCWaveFunctions/OrbitalBase.h"
#include "Optimize/VarList.h"
#include "OhmmsData/libxmldefs.h"
#include "OhmmsPETE/OhmmsVector.h"
#include "OhmmsPETE/OhmmsMatrix.h"
#include "LongRange/LRHandlerBase.h"

namespace qmcplusplus
{

class LRTwoBodyJastrow: public OrbitalBase
{
  bool NeedToRestore;
  IndexType NumPtcls;
  IndexType NumSpecies;
  IndexType NumKpts;
  ///maximum K-point
  IndexType MaxK;
  ///number of k-shells
  IndexType MaxKshell;
  /// Cell Volume
  RealType CellVolume;
  ///Rs
  RealType Rs;
  ///Normalization Constant
  RealType NormConstant;
  RealType curVal, curLap;
  PosType curGrad;

  RealType *FirstAddressOfdU;
  RealType *LastAddressOfdU;
  StructFact* skRef;
  // handler used to do evalFk
  //typedef RPALRHandler<RealType>::HandlerType HandlerType;
  typedef LRHandlerBase HandlerType;
  HandlerType* handler;

  Vector<RealType> U,d2U;
  Vector<PosType> dU;
  Vector<RealType> offU, offd2U;
  Vector<PosType> offdU;

  Matrix<ComplexType> rokbyF;
  Vector<ComplexType> Rhok;

  ////Matrix<ComplexType> rhok;
  Matrix<ComplexType> eikr;
  Vector<ComplexType> eikr_new;
  Vector<ComplexType> delta_eikr;

  vector<int> Kshell;
  vector<RealType> FkbyKK;
  //vector<PosType>  Kcart;

  void resetInternals();
  void resize();
  ///fixed components of Fk
  Vector<RealType> Fk_0;
  ///variable components of Fk
  Vector<RealType> Fk_1;
public:
  ///Fk[kindex]
  Vector<RealType> Fk;
  ///A unique Fk sorted by |k|
  Vector<RealType> Fk_symm;

  LRTwoBodyJastrow(ParticleSet& p);

  /** check out optimizable variables
   */
  void checkOutVariables(const opt_variables_type& o);

  /** check in an optimizable parameter
   * @param o a super set of optimizable variables
   */
  void checkInVariables(opt_variables_type& o);

  /** print the state, e.g., optimizables */
  void reportStatus(ostream& os);

  void resetParameters(const opt_variables_type& active);

  void resetByHandler(HandlerType* handler);

  //evaluate the distance table with els
  void resetTargetParticleSet(ParticleSet& P);

  RealType evaluateLog(ParticleSet& P,
                       ParticleSet::ParticleGradient_t& G,
                       ParticleSet::ParticleLaplacian_t& L);
  /** reset the coefficients by a function
  */
  template<typename FUNC>
  void resetByFunction(RealType kc)
  {
    RealType kcsq=kc*kc;
    int maxshell=skRef->KLists.kshell.size()-1;
    const KContainer::SContainer_t& kk(skRef->KLists.ksq);
    int ksh=0,ik=0;
    while(ksh<maxshell)
    {
      if(kk[ik]>kcsq)
        break; //exit
      ik=skRef->KLists.kshell[++ksh];
    }
    MaxKshell=ksh;
    Fk_symm.resize(MaxKshell);
    FkbyKK.resize(MaxKshell);
    Fk.resize(ik);
    Fk_0.resize(ik);
    //create a function
    FUNC uk(Rs);
    //- sign is for the form of two-body Jastrow
    RealType u0 = -4.0*M_PI/CellVolume;
    for(ksh=0,ik=0; ksh<MaxKshell; ksh++, ik++)
    {
      RealType v=u0*uk(kk[ik]);//rpa=u0/kk[ik];
      Fk_symm[ksh]=v;
      FkbyKK[ksh]=kk[ik]*v;
      for(; ik<skRef->KLists.kshell[ksh+1]; ik++)
        Fk[ik]=v;
    }
    Fk_0=Fk;
  }

  inline ValueType evaluate(ParticleSet& P,
                            ParticleSet::ParticleGradient_t& G,
                            ParticleSet::ParticleLaplacian_t& L)
  {
    return std::exp(evaluateLog(P,G,L));
  }


  ValueType ratio(ParticleSet& P, int iat);

  ValueType ratio(ParticleSet& P, int iat,
                  ParticleSet::ParticleGradient_t& dG,
                  ParticleSet::ParticleLaplacian_t& dL)
  {
    return std::exp(logRatio(P,iat,dG,dL));
  }

  ValueType logRatio(ParticleSet& P, int iat,
                     ParticleSet::ParticleGradient_t& dG,
                     ParticleSet::ParticleLaplacian_t& dL);

  void restore(int iat);
  void acceptMove(ParticleSet& P, int iat);
  void update(ParticleSet& P,
              ParticleSet::ParticleGradient_t& dG,
              ParticleSet::ParticleLaplacian_t& dL,
              int iat);


  RealType registerData(ParticleSet& P, PooledData<RealType>& buf);
  RealType updateBuffer(ParticleSet& P, PooledData<RealType>& buf, bool fromscratch=false);
  void copyFromBuffer(ParticleSet& P, PooledData<RealType>& buf);
  RealType evaluateLog(ParticleSet& P, PooledData<RealType>& buf);

  ///process input file
  bool put(xmlNodePtr cur);

  OrbitalBasePtr makeClone(ParticleSet& tqp) const;
};
}
#endif
/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$
 ***************************************************************************/
