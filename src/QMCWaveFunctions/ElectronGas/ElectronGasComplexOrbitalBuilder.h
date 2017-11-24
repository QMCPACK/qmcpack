//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Jaron T. Krogel, krogeljt@ornl.gov, Oak Ridge National Laboratory
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////
    
    
#ifndef QMCPLUSPLUS_ELECTRONGAS_COMPLEXORBITALS_H
#define QMCPLUSPLUS_ELECTRONGAS_COMPLEXORBITALS_H

#include "QMCWaveFunctions/OrbitalBuilderBase.h"
#include "QMCWaveFunctions/SPOSetBase.h"
#include "QMCWaveFunctions/BasisSetBase.h"
#include "QMCWaveFunctions/ElectronGas/HEGGrid.h"


namespace qmcplusplus
{

struct EGOSet: public SPOSetBase
{

  int KptMax;
  std::vector<PosType>  K;
  std::vector<RealType> mK2;

  EGOSet(const std::vector<PosType>& k, const std::vector<RealType>& k2);
  EGOSet(const std::vector<PosType>& k, const std::vector<RealType>& k2, const std::vector<int>& d);

  SPOSetBase* makeClone() const
  {
    return new EGOSet(*this);
  }

  void resetParameters(const opt_variables_type& optVariables) {}
  inline void resetTargetParticleSet(ParticleSet& P) { }
  void setOrbitalSetSize(int norbs) { }

  inline void
  evaluate(const ParticleSet& P, int iat, ValueVector_t& psi)
  {
    const PosType &r=P.activeR(iat);
    RealType sinkr,coskr;
    for(int ik=0; ik<KptMax; ik++)
    {
      sincos(dot(K[ik],r),&sinkr,&coskr);
      psi[ik]=ValueType(coskr,sinkr);
    }
  }

  /** generic inline function to handle a row
   * @param P current ParticleSet
   * @param iat active particle
   * @param psi value row
   * @param dpsi gradient row
   * @param d2psi laplacian row
   */
  inline void
  evaluate(const ParticleSet& P, int iat,
           ValueVector_t& psi, GradVector_t& dpsi, ValueVector_t& d2psi)
  {
    const PosType &r=P.activeR(iat);
    RealType sinkr,coskr;
    for(int ik=0; ik<KptMax; ik++)
    {
      sincos(dot(K[ik],r),&sinkr,&coskr);
      psi[ik]  =ValueType(coskr,sinkr);
      dpsi[ik] =ValueType(-sinkr,coskr)*K[ik];
      d2psi[ik]=ValueType(mK2[ik]*coskr,mK2[ik]*sinkr);
    }
  }

  void evaluate(const ParticleSet& P, int iat, ValueVector_t& psi, GradVector_t& dpsi
                , HessVector_t& grad_grad_psi)
  {
    APP_ABORT("Incomplete implementation EGOSet::evaluate(P,iat,psi,dpsi,grad_grad_psi)");
  }

  void evaluate_notranspose(const ParticleSet& P, int first, int last,
                            ValueMatrix_t& logdet, GradMatrix_t& dlogdet, ValueMatrix_t& d2logdet)
  {
    for(int iat=first, i=0; iat<last; ++iat,++i)
    {
      ValueVector_t v(logdet[i],OrbitalSetSize);
      GradVector_t g(dlogdet[i],OrbitalSetSize);
      ValueVector_t l(d2logdet[i],OrbitalSetSize);
      evaluate(P,iat,v,g,l);
    }
  }
  void evaluate_notranspose(const ParticleSet& P, int first, int last
                            , ValueMatrix_t& logdet, GradMatrix_t& dlogdet, HessMatrix_t& grad_grad_logdet, GGGMatrix_t& grad_grad_grad_logdet)
  {
    APP_ABORT("Incomplete implementation EGOSet::evaluate(P,first,last,lodget,dlodet,grad_grad_logdet,grad_grad_grad_logdet");
  }
};


/** OrbitalBuilder for Slater determinants of electron-gas
*/
class ElectronGasComplexOrbitalBuilder: public OrbitalBuilderBase
{
public:

  ///constructor
  ElectronGasComplexOrbitalBuilder(ParticleSet& els, TrialWaveFunction& wfs);
  //typedef VarRegistry<RealType> OptimizableSetType;

  ///implement vritual function
  bool put(xmlNodePtr cur);

};

/** OrbitalBuilder for Slater determinants of electron-gas
*/
class ElectronGasBasisBuilder: public BasisSetBuilder
{
protected:
  bool has_twist;
  PosType unique_twist;
  HEGGrid<RealType,OHMMS_DIM> egGrid;
  xmlNodePtr spo_node;
public:
  ///constructor
  ElectronGasBasisBuilder(ParticleSet& p, xmlNodePtr cur);

  ///implement virtual function
  bool put(xmlNodePtr cur);
  /** initialize the Antisymmetric wave function for electrons
  *@param cur the current xml node
  */
  SPOSetBase* createSPOSetFromXML(xmlNodePtr cur);
  SPOSetBase* createSPOSetFromIndices(indices_t& indices);

};
}
#endif
