//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Miguel Morales, moralessilva2@llnl.gov, Lawrence Livermore National Laboratory
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Jaron T. Krogel, krogeljt@ornl.gov, Oak Ridge National Laboratory
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////
    
    
#ifndef QMCPLUSPLUS_ELECTRONGAS_ORBITALS_H
#define QMCPLUSPLUS_ELECTRONGAS_ORBITALS_H

#include <QMCWaveFunctions/OrbitalBuilderBase.h>
#include <QMCWaveFunctions/SPOSetBase.h>
#include <config/stdlib/math.h>

#include "QMCWaveFunctions/BasisSetBase.h"
#include "QMCWaveFunctions/ElectronGas/HEGGrid.h"

#if defined(QMC_COMPLEX)
#include "QMCWaveFunctions/ElectronGas/ElectronGasComplexOrbitalBuilder.h"
#else /** declare real HEG orbitals **/
namespace qmcplusplus
{

//forward declaration
class  BackflowTransformation;

struct RealEGOSet: public SPOSetBase
{

  int KptMax;
  RealType kdotr;
  std::vector<PosType> K;
  std::vector<RealType> mK2;

  RealEGOSet(const std::vector<PosType>& k, const std::vector<RealType>& k2);

  void resetParameters(const opt_variables_type& optVariables) {}
  inline void resetTargetParticleSet(ParticleSet& P) { }
  void setOrbitalSetSize(int norbs) { }

  SPOSetBase* makeClone() const
  {
    return new RealEGOSet(*this);
  }

  PosType get_k(int i)
  {
    //Only used in the same_k part of the optimizable SPO set. we allow optimization to k points in the same direction
    if(i>0)
    {
      int ik=(i-1)/2;
      int even=(i-1)%2;
      PosType k_tmp = K[ik];
      k_tmp *= 1.0/std::sqrt(-mK2[ik]);
//         if(even)
//           return -1*k_tmp;
//         else
      return k_tmp;
    }
    else
      return PosType();
  };

  inline ValueType f(const PosType& pos,int i)
  {
    if(i>0)
    {
      int ik=(i-1)/2;
      int even=(i-1)%2;
      kdotr=dot(K[ik],pos);
      if(even)
        return std::cos(kdotr);
      else
        return std::sin(kdotr);
    }
    else
      return 1.0;
  }

  void evaluate(const ParticleSet& P, int iat, ValueVector_t& psi)
  {
    const PosType &r=P.activeR(iat);
    RealType sinkr,coskr;
    psi[0]=1.0;
    for(int ik=0, j=1; ik<KptMax; ik++)
    {
      sincos(dot(K[ik],r),&sinkr,&coskr);
      psi[j++]=coskr;
      psi[j++]=sinkr;
    }
  }

  /** generic inline function to handle a row
   * @param P current ParticleSet
   * @param iat active particle
   * @param psi value row
   * @param dpsi gradient row
   * @param d2psi laplacian row
   */
  inline void evaluate(const ParticleSet& P, int iat, ValueVector_t& psi, GradVector_t& dpsi, ValueVector_t& d2psi)
  {
    psi[0]=1.0;
    dpsi[0]=0.0;
    d2psi[0]=0.0;
    const PosType &r=P.activeR(iat);
    RealType coskr, sinkr;
    for(int ik=0,j1=1; ik<KptMax; ik++,j1+=2)
    {
      int j2=j1+1;
      sincos(dot(K[ik],r),&sinkr,&coskr);
      psi[j1]=coskr;
      psi[j2]=sinkr;
      dpsi[j1]=-sinkr*K[ik];
      dpsi[j2]= coskr*K[ik];
      d2psi[j1]=mK2[ik]*coskr;
      d2psi[j2]=mK2[ik]*sinkr;
    }
  }

  /** generic inline function to handle a row
   * @param P current ParticleSet
   * @param iat active particle
   * @param psi value row
   * @param dpsi gradient row
   * @param hess hessian row
   */
  inline void evaluate(const ParticleSet& P, int iat, ValueVector_t& psi, GradVector_t& dpsi, HessVector_t& hess)
  {
    psi[0]=1.0;
    dpsi[0]=0.0;
    hess[0]=0.0;
    const PosType &r=P.activeR(iat);
    RealType coskr, sinkr;
    for(int ik=0,j1=1; ik<KptMax; ik++,j1+=2)
    {
      int j2=j1+1;
      sincos(dot(K[ik],r),&sinkr,&coskr);
      psi[j1]=coskr;
      psi[j2]=sinkr;
      dpsi[j1]=-sinkr*K[ik];
      dpsi[j2]= coskr*K[ik];
      for(int la=0; la<OHMMS_DIM; la++)
      {
        (hess[j1])(la,la)=-coskr*(K[ik])[la]*(K[ik])[la];
        (hess[j2])(la,la)=-sinkr*(K[ik])[la]*(K[ik])[la];
        for(int lb=la+1; lb<OHMMS_DIM; lb++)
        {
          (hess[j1])(la,lb)=-coskr*(K[ik])[la]*(K[ik])[lb];
          (hess[j2])(la,lb)=-sinkr*(K[ik])[la]*(K[ik])[lb];
          (hess[j1])(lb,la)=(hess[j1])(la,lb);
          (hess[j2])(lb,la)=(hess[j2])(la,lb);
        }
      }
    }
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
                            , ValueMatrix_t& logdet, GradMatrix_t& dlogdet
                            , HessMatrix_t& grad_grad_logdet)
  {
    for(int i=0,iat=first; iat<last; i++,iat++)
    {
      ValueType* psi=logdet[i];
      GradType* dpsi=dlogdet[i];
      HessType*  hess=grad_grad_logdet[i];
      psi[0]=1.0;
      dpsi[0]=0.0;
      hess[0]=0.0;
      RealType coskr, sinkr;
      for(int ik=0,j1=1; ik<KptMax; ik++,j1+=2)
      {
        int j2=j1+1;
        sincos(dot(K[ik],P.R[iat]),&sinkr,&coskr);
        psi[j1]=coskr;
        psi[j2]=sinkr;
        dpsi[j1]=-sinkr*K[ik];
        dpsi[j2]= coskr*K[ik];
        for(int la=0; la<OHMMS_DIM; la++)
        {
          (hess[j1])(la,la)=-coskr*(K[ik])[la]*(K[ik])[la];
          (hess[j2])(la,la)=-sinkr*(K[ik])[la]*(K[ik])[la];
          for(int lb=la+1; lb<OHMMS_DIM; lb++)
          {
            (hess[j1])(la,lb)=-coskr*(K[ik])[la]*(K[ik])[lb];
            (hess[j2])(la,lb)=-sinkr*(K[ik])[la]*(K[ik])[lb];
            (hess[j1])(lb,la)=(hess[j1])(la,lb);
            (hess[j2])(lb,la)=(hess[j2])(la,lb);
          }
        }
      }
    }
  }

  void evaluate_notranspose(const ParticleSet& P, int first, int last
                            , ValueMatrix_t& logdet, GradMatrix_t& dlogdet, HessMatrix_t& grad_grad_logdet, GGGMatrix_t& grad_grad_grad_logdet)
  {
    for(int i=0,iat=first; iat<last; i++,iat++)
    {
      ValueType* psi=logdet[i];
      GradType* dpsi=dlogdet[i];
      HessType*  hess=grad_grad_logdet[i];
      GGGType* ggg=grad_grad_grad_logdet[i];
      psi[0]=1.0;
      dpsi[0]=0.0;
      hess[0]=0.0;
      ggg[0]=0.0;
      RealType coskr, sinkr;
      for(int ik=0,j1=1; ik<KptMax; ik++,j1+=2)
      {
        int j2=j1+1;
        sincos(dot(K[ik],P.R[iat]),&sinkr,&coskr);
        psi[j1]=coskr;
        psi[j2]=sinkr;
        dpsi[j1]=-sinkr*K[ik];
        dpsi[j2]= coskr*K[ik];
        for(int la=0; la<OHMMS_DIM; la++)
        {
          (hess[j1])(la,la)=-coskr*(K[ik])[la]*(K[ik])[la];
          (hess[j2])(la,la)=-sinkr*(K[ik])[la]*(K[ik])[la];
          ((ggg[j1])[la])(la,la)=sinkr*(K[ik])[la]*(K[ik])[la]*(K[ik])[la];
          ((ggg[j2])[la])(la,la)=-coskr*(K[ik])[la]*(K[ik])[la]*(K[ik])[la];
          for(int lb=la+1; lb<OHMMS_DIM; lb++)
          {
            (hess[j1])(la,lb)=-coskr*(K[ik])[la]*(K[ik])[lb];
            (hess[j2])(la,lb)=-sinkr*(K[ik])[la]*(K[ik])[lb];
            (hess[j1])(lb,la)=(hess[j1])(la,lb);
            (hess[j2])(lb,la)=(hess[j2])(la,lb);
            ((ggg[j1])[la])(lb,la) = sinkr*(K[ik])[la]*(K[ik])[lb]*(K[ik])[la];
            ((ggg[j2])[la])(lb,la) = -coskr*(K[ik])[la]*(K[ik])[lb]*(K[ik])[la];
            ((ggg[j1])[la])(la,lb) = ((ggg[j1])[la])(lb,la);
            ((ggg[j2])[la])(la,lb) = ((ggg[j2])[la])(lb,la);
            ((ggg[j1])[lb])(la,la) = ((ggg[j1])[la])(lb,la);
            ((ggg[j2])[lb])(la,la) = ((ggg[j2])[la])(lb,la);
            ((ggg[j1])[la])(lb,lb) = sinkr*(K[ik])[la]*(K[ik])[lb]*(K[ik])[lb];
            ((ggg[j2])[la])(lb,lb) = -coskr*(K[ik])[la]*(K[ik])[lb]*(K[ik])[lb];
            ((ggg[j1])[lb])(la,lb) = ((ggg[j1])[la])(lb,lb);
            ((ggg[j2])[lb])(la,lb) = ((ggg[j2])[la])(lb,lb);
            ((ggg[j1])[lb])(lb,la) = ((ggg[j1])[la])(lb,lb);
            ((ggg[j2])[lb])(lb,la) = ((ggg[j2])[la])(lb,lb);
            for(int lc=lb+1; lc<OHMMS_DIM; lc++)
            {
              ( (ggg[j1])[la] )(lb,lc) = sinkr*(K[ik])[la]*(K[ik])[lb]*(K[ik])[lc];
              ( (ggg[j2])[la] )(lb,lc) = -coskr*(K[ik])[la]*(K[ik])[lb]*(K[ik])[lc];
              ( (ggg[j1])[la] )(lc,lb) = ( (ggg[j1])[la] )(lb,lc);
              ( (ggg[j2])[la] )(lc,lb) = ( (ggg[j2])[la] )(lb,lc);
              ( (ggg[j1])[lb] )(la,lc) = ( (ggg[j1])[la] )(lb,lc);
              ( (ggg[j2])[lb] )(la,lc) = ( (ggg[j2])[la] )(lb,lc);
              ( (ggg[j1])[lb] )(lc,la) = ( (ggg[j1])[la] )(lb,lc);
              ( (ggg[j2])[lb] )(lc,la) = ( (ggg[j2])[la] )(lb,lc);
              ( (ggg[j1])[lc] )(la,lb) = ( (ggg[j1])[la] )(lb,lc);
              ( (ggg[j2])[lc] )(la,lb) = ( (ggg[j2])[la] )(lb,lc);
              ( (ggg[j1])[lc] )(lb,la) = ( (ggg[j1])[la] )(lb,lc);
              ( (ggg[j2])[lc] )(lb,la) = ( (ggg[j2])[la] )(lb,lc);
            }
          }
        }
//#if OHMMS_DIM ==3
//          for(int la=0; la<3; la++) {
//            for(int lb=0; lb<3; lb++) {
//              for(int lc=0; lc<3; lc++) {
//                ( (ggg[j1])[la] )(lb,lc) = sinkr*(K[ik])[la]*(K[ik])[lb]*(K[ik])[lc];
//                ( (ggg[j2])[la] )(lb,lc) = -coskr*(K[ik])[la]*(K[ik])[lb]*(K[ik])[lc];
//              }
//            }
//          }
//#elif OHMMS_DIM ==2
//          for(int la=0; la<2; la++) {
//            for(int lb=0; lb<2; lb++) {
//              for(int lc=0; lc<2; lc++) {
//                ( (ggg[j1])[la] )(lb,lc) = sinkr*(K[ik])[la]*(K[ik])[lb]*(K[ik])[lc];
//                ( (ggg[j2])[la] )(lb,lc) = -coskr*(K[ik])[la]*(K[ik])[lb]*(K[ik])[lc];
//              }
//            }
//          }
//#endif
      }
    }
  }


};

/** OrbitalBuilder for Slater determinants of electron-gas
*/
class ElectronGasOrbitalBuilder: public OrbitalBuilderBase
{
public:

  ///constructor
  ElectronGasOrbitalBuilder(ParticleSet& els, TrialWaveFunction& wfs);

  ///implement vritual function
  bool put(xmlNodePtr cur);

  bool UseBackflow;
  BackflowTransformation *BFTrans;
};

/** OrbitalBuilder for Slater determinants of electron-gas
*/
class ElectronGasBasisBuilder: public BasisSetBuilder
{
protected:
  HEGGrid<RealType,OHMMS_DIM> egGrid;
  xmlNodePtr spo_node;
public:
  ///constructor
  ElectronGasBasisBuilder(ParticleSet& p, xmlNodePtr cur);

  ///implement vritual function
  bool put(xmlNodePtr cur);
  /** initialize the Antisymmetric wave function for electrons
  *@param cur the current xml node
  */
  SPOSetBase* createSPOSetFromXML(xmlNodePtr cur);
};

}
#endif /** QMC_COMPLEX **/
#endif
