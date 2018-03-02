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
    
    
#ifndef QMCPLUSPLUS_BACKFLOW_ELEC_ION_SPINH
#define QMCPLUSPLUS_BACKFLOW_ELEC_ION_SPINH
#include "QMCWaveFunctions/OrbitalSetTraits.h"
#include "QMCWaveFunctions/Fermion/BackflowFunctionBase.h"
#include "Particle/DistanceTable.h"
#include <cmath>
#include <vector>

namespace qmcplusplus
{
template<class FT>
class Backflow_eI_spin: public BackflowFunctionBase
{

public:

  bool Spin;
  Matrix<FT*> RadFunc;
  Matrix<int> Fmask;
  Matrix<int> offsetPrms;
  /** index offset for the source particles
   *
   * The loop over the particles is replaced by a double loop as
   * for(int sg=0;sg<F.rows();++sg)
   *   for(int s=s_offset[sg]; s<s_offset[sg+1];++s)
   */
  std::vector<int> s_offset;
  /** index offset for the target particles
   *
   * The loop over the particles is replaced by a double loop as
   * for(int tg=0;tg<F.cols();++tg)
   *   for(int t=t_offset[tg]; t<t_offset[tg+1];++t)
   */
  std::vector<int> t_offset;

  Backflow_eI_spin(ParticleSet& ions, ParticleSet& els)
    : BackflowFunctionBase(ions,els), Spin(false)
  {
    myTable = DistanceTable::add(ions,els,DT_AOS);
    RadFunc.resize(ions.groups(),els.groups());
    for(int i=0; i<RadFunc.size(); ++i)
      RadFunc(i)=0;
    Fmask.resize(ions.groups(),els.groups());
    Fmask=-1;
    resize(NumTargets,NumCenters);
    s_offset.resize(ions.groups()+1,0);
    t_offset.resize(els.groups()+1,0);
    for(int s=0; s<ions.groups(); ++s)
      s_offset[s+1]=ions.last(s);
    for(int t=0; t<els.groups(); ++t)
      t_offset[t+1]=els.last(t);
    offsetPrms.resize(ions.groups(),els.groups());
    offsetPrms=0;
  }

  ~Backflow_eI_spin()
  {
    for(int i=0; i<RadFunc.size(); ++i)
      if(RadFunc[i])
        delete RadFunc[i];
  }

  void addFunc(int source_g, FT* afunc, int target_g)
  {
    if(target_g<0)
    {
      if(Spin)
      {
        APP_ABORT("Cannot mix spin-dependent with spin-indepdentent Jastrow");
      }
      int pid=source_g*RadFunc.cols();
      for(int ig=0; ig<RadFunc.cols(); ++ig)
      {
        RadFunc(source_g,ig)=afunc;
        Fmask(source_g,ig)=pid;
      }
      app_log() << " Adding functor of type "  << source_g << " for all the target. " << std::endl;
    }
    else
    {
      Spin=true;
      RadFunc(source_g,target_g)=afunc;
      Fmask(source_g,target_g)=source_g*RadFunc.cols()+target_g;
      app_log() << " Adding functor of type "  << source_g << " for the target type " << target_g << std::endl;
    }
  }

  void resetTargetParticleSet(ParticleSet& P)
  {
    myTable = DistanceTable::add(CenterSys,P,DT_AOS);
  }

  BackflowFunctionBase* makeClone(ParticleSet& tqp)
  {
    Backflow_eI_spin<FT>* clone = new Backflow_eI_spin<FT>(CenterSys,tqp);
    clone->resize(NumTargets,NumCenters);
    clone->indexOfFirstParam=indexOfFirstParam;
    clone->offsetPrms=offsetPrms;
    clone->numParams=numParams;
    clone->derivs=derivs;
    if(Spin)
    {
      for(int sg=0; sg<RadFunc.rows(); ++sg)
        for(int tg=0; tg<RadFunc.rows(); ++tg)
          if(RadFunc(sg,tg))
            clone->addFunc(sg,new FT(*RadFunc(sg,tg)),tg);
    }
    else
    {
      for(int sg=0; sg<RadFunc.rows(); ++sg)
        if(RadFunc(sg,0))
          clone->addFunc(sg,new FT(*RadFunc(sg,0)),-1);
    }
    return clone;
  }

  void reportStatus(std::ostream& os)
  {
    std::cerr <<RadFunc.size()<< std::endl;
    std::cerr <<isOptimizable()<< std::endl;
    std::cerr <<myTable<< std::endl;
    std::cerr << offsetPrms << std::endl;
  }

  void resetParameters(const opt_variables_type& active)
  {
    for(int i=0; i<RadFunc.size(); ++i)
      if(Fmask(i) == i)
        RadFunc(i)->resetParameters(active);
  }

  void checkInVariables(opt_variables_type& active)
  {
    //myVars.clear();
    for(int i=0; i<RadFunc.size(); ++i)
      if(Fmask(i) == i)
        RadFunc(i)->checkInVariables(active);
  }

  void checkOutVariables(const opt_variables_type& active)
  {
    for(int i=0; i<RadFunc.size(); ++i)
      if(Fmask(i) == i)
        RadFunc(i)->checkOutVariables(active);
  }

  inline bool isOptimizable()
  {
    for(int i=0; i<RadFunc.size(); ++i)
      if(Fmask(i) == i && RadFunc(i)->isOptimizable())
        return true;
    return false;
  }

  inline int
  indexOffset()
  {
    for(int i=0; i<RadFunc.size(); ++i)
      if(Fmask(i) == i)
        return RadFunc(i)->myVars.where(0);
    return -1;
  }

  // only elements of qp coordinates that changed should be copied
  inline void
  acceptMove(int iat, int UpdateMode)
  {
    int num;
    switch(UpdateMode)
    {
    case ORB_PBYP_RATIO:
      num = UIJ.cols();
      for(int i=0; i<num; i++)
        UIJ(iat,i) = UIJ_temp[i];
      break;
    case ORB_PBYP_PARTIAL:
      num = UIJ.cols();
      for(int i=0; i<num; i++)
        UIJ(iat,i) = UIJ_temp[i];
      num = AIJ.cols();
      for(int i=0; i<num; i++)
        AIJ(iat,i) = AIJ_temp[i];
      break;
    case ORB_PBYP_ALL:
      num = UIJ.cols();
      for(int i=0; i<num; i++)
        UIJ(iat,i) = UIJ_temp[i];
      num = AIJ.cols();
      for(int i=0; i<num; i++)
        AIJ(iat,i) = AIJ_temp[i];
      num = BIJ.cols();
      for(int i=0; i<num; i++)
        BIJ(iat,i) = BIJ_temp[i];
      break;
    default:
      num = UIJ.cols();
      for(int i=0; i<num; i++)
        UIJ(iat,i) = UIJ_temp[i];
      num = AIJ.cols();
      for(int i=0; i<num; i++)
        AIJ(iat,i) = AIJ_temp[i];
      num = BIJ.cols();
      for(int i=0; i<num; i++)
        BIJ(iat,i) = BIJ_temp[i];
      break;
    }
    UIJ_temp=0.0;
    AIJ_temp=0.0;
    BIJ_temp=0.0;
  }

  inline void
  restore(int iat, int UpdateType)
  {
    UIJ_temp=0.0;
    AIJ_temp=0.0;
    BIJ_temp=0.0;
  }

  void registerData(WFBufferType& buf)
  {
    FirstOfU = &(UIJ(0,0)[0]);
    LastOfU = FirstOfU + OHMMS_DIM*NumTargets*NumCenters;
    FirstOfA = &(AIJ(0,0)[0]);
    LastOfA = FirstOfA + OHMMS_DIM*OHMMS_DIM*NumTargets*NumCenters;
    FirstOfB = &(BIJ(0,0)[0]);
    LastOfB = FirstOfB + OHMMS_DIM*NumTargets*NumCenters;
    buf.add(FirstOfU,LastOfU);
    buf.add(FirstOfA,LastOfA);
    buf.add(FirstOfB,LastOfB);
  }

  /** calculate quasi-particle coordinates only
   */
  inline void
  evaluate(const ParticleSet& P, ParticleSet& QP)
  {
    RealType du,d2u;
    for(int sg=0; sg<RadFunc.rows(); ++sg)
    {
      for(int iat=s_offset[sg]; iat< s_offset[sg+1]; ++iat)
      {
        int nn=myTable->M[iat];//starting nn for the iat-th source
        for(int tg=0; tg<RadFunc.cols(); ++tg)
        {
          FT* func=RadFunc(sg,tg);
          if(func)
            for(int jat=t_offset[tg]; jat< t_offset[tg+1]; ++jat,++nn)
            {
              RealType uij = func->evaluate(myTable->r(nn),du,d2u);
              QP.R[jat] += (UIJ(jat,iat) = uij*myTable->dr(nn));  // dr(ij) = r_j-r_i
            }
          else
            nn+=t_offset[tg+1]-t_offset[tg];//move forward by the number of particles in the group tg
        }
      }
    }
  }


  inline void
  evaluate(const ParticleSet& P, ParticleSet& QP, GradVector_t& Bmat, HessMatrix_t& Amat)
  {
    RealType du,d2u,temp;
    for(int sg=0; sg<RadFunc.rows(); ++sg)
    {
      for(int iat=s_offset[sg]; iat< s_offset[sg+1]; ++iat)
      {
        int nn=myTable->M[iat];//starting nn for the iat-th source
        for(int tg=0; tg<RadFunc.cols(); ++tg)
        {
          FT* func=RadFunc(sg,tg);
          if(func)
            for(int jat=t_offset[tg]; jat< t_offset[tg+1]; ++jat,++nn)
            {
              RealType uij = func->evaluate(myTable->r(nn),du,d2u);
              //PosType u = uij*myTable->dr(nn);
              du *= myTable->rinv(nn);
              QP.R[jat] += (UIJ(jat,iat) = uij*myTable->dr(nn));
              HessType& hess = AIJ(jat,iat);
              hess = du*outerProduct(myTable->dr(nn),myTable->dr(nn));
              hess[0] += uij;
              hess[4] += uij;
              hess[8] += uij;
              Amat(jat,jat) += hess;
              //u = (d2u+4.0*du)*myTable->dr(nn);
              Bmat[jat] += (BIJ(jat,iat)=(d2u+4.0*du)*myTable->dr(nn));
            }
          else
            nn+=t_offset[tg+1]-t_offset[tg];//move forward by the number of particles in the group tg
        }
      }
    }
  }


  /** calculate quasi-particle coordinates, Bmat and Amat
   */
  inline void
  evaluate(const ParticleSet& P, ParticleSet& QP, GradMatrix_t& Bmat_full, HessMatrix_t& Amat)
  {
    RealType du,d2u,temp;
    for(int sg=0; sg<RadFunc.rows(); ++sg)
    {
      for(int iat=s_offset[sg]; iat< s_offset[sg+1]; ++iat)
      {
        int nn=myTable->M[iat];//starting nn for the iat-th source
        for(int tg=0; tg<RadFunc.cols(); ++tg)
        {
          FT* func=RadFunc(sg,tg);
          if(func)
            for(int jat=t_offset[tg]; jat< t_offset[tg+1]; ++jat,++nn)
            {
              RealType uij = func->evaluate(myTable->r(nn),du,d2u);
              du *= myTable->rinv(nn);
              //PosType u = uij*myTable->dr(nn);
              QP.R[jat] += (UIJ(jat,iat) = uij*myTable->dr(nn));
              HessType& hess = AIJ(jat,iat);
              hess = du*outerProduct(myTable->dr(nn),myTable->dr(nn));
              hess[0] += uij;
              hess[4] += uij;
              hess[8] += uij;
              Amat(jat,jat) += hess;
              // this will create problems with QMC_COMPLEX, because Bmat is ValueType and dr is RealType
              //u = (d2u+4.0*du)*myTable->dr(nn);
              Bmat_full(jat,jat) += (BIJ(jat,iat)=(d2u+4.0*du)*myTable->dr(nn));
            }
          else
            nn+=t_offset[tg+1]-t_offset[tg];//move forward by the number of particles in the group tg
        }
      }
    }
  }

  /** calculate quasi-particle coordinates after pbyp move
   */
  inline void
  evaluatePbyP(const ParticleSet& P, ParticleSet::ParticlePos_t& newQP
               ,const std::vector<int>& index)
  {
    evaluatePbyP(P,index[0],newQP);
  }


  /** calculate quasi-particle coordinates after pbyp move
   */
  inline void
  evaluatePbyP(const ParticleSet& P, int iat, ParticleSet::ParticlePos_t& newQP)
  {
    RealType du,d2u;
    int tg=P.GroupID[iat];//species of this particle
    for(int sg=0; sg<RadFunc.rows(); ++sg)
    {
      FT* func=RadFunc(sg,tg);
      if(func)
      {
        for(int j=s_offset[sg]; j< s_offset[sg+1]; ++j)
        {
          RealType uij = func->evaluate(myTable->Temp[j].r1,du,d2u);
          PosType u = (UIJ_temp[j]=uij*myTable->Temp[j].dr1)-UIJ(iat,j);
          newQP[iat] += u;
        }
      }
    }
  }

  inline void
  evaluatePbyP(const ParticleSet& P, ParticleSet::ParticlePos_t& newQP
               ,const std::vector<int>& index, HessMatrix_t& Amat)
  {
    evaluatePbyP(P,index[0],newQP,Amat);
  }

  inline void
  evaluatePbyP(const ParticleSet& P, int iat
               ,ParticleSet::ParticlePos_t& newQP, HessMatrix_t& Amat)
  {
    RealType du,d2u;
    int tg=P.GroupID[iat];//species of this particle
    for(int sg=0; sg<RadFunc.rows(); ++sg)
    {
      FT* func=RadFunc(sg,tg);
      if(func)
      {
        for(int j=s_offset[sg]; j< s_offset[sg+1]; ++j)
        {
          RealType uij = func->evaluate(myTable->Temp[j].r1,du,d2u);
          PosType u = (UIJ_temp[j]=uij*myTable->Temp[j].dr1)-UIJ(iat,j);
          newQP[iat] += u;
          HessType& hess = AIJ_temp[j];
          hess = (du*myTable->Temp[j].rinv1)*outerProduct(myTable->Temp[j].dr1,myTable->Temp[j].dr1);
          hess[0] += uij;
          hess[4] += uij;
          hess[8] += uij;
          // should I expand this??? Is the compiler smart enough???
          Amat(iat,iat) += (hess - AIJ(iat,j));
        }
      }
    }
  }

  inline void
  evaluatePbyP(const ParticleSet& P, ParticleSet::ParticlePos_t& newQP
               ,const std::vector<int>& index, GradMatrix_t& Bmat_full, HessMatrix_t& Amat)
  {
    evaluatePbyP(P,index[0],newQP,Bmat_full,Amat);
  }

  inline void
  evaluatePbyP(const ParticleSet& P, int iat, ParticleSet::ParticlePos_t& newQP
               , GradMatrix_t& Bmat_full, HessMatrix_t& Amat)
  {
    RealType du,d2u;
    int tg=P.GroupID[iat];//species of this particle
    for(int sg=0; sg<RadFunc.rows(); ++sg)
    {
      FT* func=RadFunc(sg,tg);
      if(func)
      {
        for(int j=s_offset[sg]; j< s_offset[sg+1]; ++j)
        {
          RealType uij = func->evaluate(myTable->Temp[j].r1,du,d2u);
          PosType u = (UIJ_temp[j]=uij*myTable->Temp[j].dr1)-UIJ(iat,j);
          newQP[iat] += u;
          du *= myTable->Temp[j].rinv1;
          HessType& hess = AIJ_temp[j];
          hess = du*outerProduct(myTable->Temp[j].dr1,myTable->Temp[j].dr1);
          hess[0] += uij;
          hess[4] += uij;
          hess[8] += uij;
          Amat(iat,iat) += (hess - AIJ(iat,j));
          BIJ_temp[j]=(d2u+4.0*du)*myTable->Temp[j].dr1;
          Bmat_full(iat,iat) += (BIJ_temp[j]-BIJ(iat,j));
        }
      }
    }
  }

  /** calculate only Bmat
   *  This is used in pbyp moves, in updateBuffer()
   */
  inline void
  evaluateBmatOnly(const ParticleSet& P,GradMatrix_t& Bmat_full)
  {
    RealType du,d2u;
    for(int sg=0; sg<RadFunc.rows(); ++sg)
    {
      for(int iat=s_offset[sg]; iat< s_offset[sg+1]; ++iat)
      {
        int nn=myTable->M[iat];//starting nn for the iat-th source
        for(int tg=0; tg<RadFunc.cols(); ++tg)
        {
          FT* func=RadFunc(sg,tg);
          if(func)
            for(int jat=t_offset[tg]; jat< t_offset[tg+1]; ++jat,++nn)
            {
              RealType uij = func->evaluate(myTable->r(nn),du,d2u);
              Bmat_full(jat,jat) += (BIJ(jat,iat)=(d2u+4.0*du*myTable->rinv(nn))*myTable->dr(nn));
            }
          else
            nn+=t_offset[tg+1]-t_offset[tg];//move forward by the number of particles in the group tg
        }
      }
    }
  }

  /** calculate quasi-particle coordinates, Bmat and Amat
   *  calculate derivatives wrt to variational parameters
   */
  inline void
  evaluateWithDerivatives(const ParticleSet& P, ParticleSet& QP, GradMatrix_t& Bmat_full, HessMatrix_t& Amat, GradMatrix_t& Cmat, GradMatrix_t& Ymat, HessArray_t& Xmat)
  {
    RealType du,d2u,temp;
    for(int sg=0; sg<RadFunc.rows(); ++sg)
    {
      for(int iat=s_offset[sg]; iat< s_offset[sg+1]; ++iat)
      {
        int nn=myTable->M[iat];//starting nn for the iat-th source
        for(int tg=0; tg<RadFunc.cols(); ++tg)
        {
          FT* func=RadFunc(sg,tg);
          if(func)
            for(int jat=t_offset[tg]; jat< t_offset[tg+1]; ++jat,++nn)
            {
              RealType uij = func->evaluate(myTable->r(nn),du,d2u);
              //           std::fill(derivs.begin(),derivs.end(),0.0);
              int NPrms = func->NumParams;
              std::vector<TinyVector<RealType,3> > derivsju(NPrms);
              func->evaluateDerivatives(myTable->r(nn),derivsju);
              du *= myTable->rinv(nn);
              //PosType u = uij*myTable->dr(nn);
              QP.R[jat] += (UIJ(jat,iat) = uij*myTable->dr(nn));
              HessType op = outerProduct(myTable->dr(nn),myTable->dr(nn));
              HessType& hess = AIJ(jat,iat);
              hess = du*op;
              hess[0] += uij;
              hess[4] += uij;
              hess[8] += uij;
              Amat(jat,jat) += hess;
              // this will create problems with QMC_COMPLEX, because Bmat is ValueType and dr is RealType
              //u = (d2u+4.0*du)*myTable->dr(nn);
              Bmat_full(jat,jat) += (BIJ(jat,iat)=(d2u+4.0*du)*myTable->dr(nn));
              //WARNINGL offsetPrms
              for(int prm=0,la=indexOfFirstParam+offsetPrms(sg,tg); prm<NPrms; prm++,la++)
              {
                Cmat(la,jat) += myTable->dr(nn)*derivsju[prm][0];
                Xmat(la,jat,jat) += (derivsju[prm][1]*myTable->rinv(nn))*op;
                Xmat(la,jat,jat)[0] += derivsju[prm][0];
                Xmat(la,jat,jat)[4] += derivsju[prm][0];
                Xmat(la,jat,jat)[8] += derivsju[prm][0];
                Ymat(la,jat) += (derivsju[prm][2]+4.0*derivsju[prm][1]*myTable->rinv(nn))*myTable->dr(nn);
              }
            }
          else
            nn+=t_offset[tg+1]-t_offset[tg];//move forward by the number of particles in the group tg
        }
      }
    }
  }

};
}

#endif
