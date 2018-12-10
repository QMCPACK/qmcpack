//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Miguel Morales, moralessilva2@llnl.gov, Lawrence Livermore National Laboratory
//                    Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Miguel Morales, moralessilva2@llnl.gov, Lawrence Livermore National Laboratory
//////////////////////////////////////////////////////////////////////////////////////
    
    
#ifndef QMCPLUSPLUS_BACKFLOW_EE_KSPACE_H
#define QMCPLUSPLUS_BACKFLOW_EE_KSPACE_H
#include "QMCWaveFunctions/OrbitalSetTraits.h"
#include "QMCWaveFunctions/Fermion/BackflowFunctionBase.h"
#include "Particle/DistanceTable.h"
#include <LongRange/StructFact.h>
#include "Message/Communicate.h"
#include <cmath>

namespace qmcplusplus
{
class Backflow_ee_kSpace: public BackflowFunctionBase
{

  typedef QMCTraits::ComplexType   ComplexType;
  ///typedef for real values
  //typedef optimize::VariableSet::real_type real_type;
  ///typedef for variableset: this is going to be replaced
  typedef optimize::VariableSet opt_variables_type;
  ///typedef for name-value lists
  typedef optimize::VariableSet::variable_map_type variable_map_type;

public:

  //number of groups of the target particleset
  bool Optimize;
  int numParams;
  std::vector<RealType> Fk;
  std::vector<int> offsetPrms;
  int NumGroups;
  int NumKShells; // number of k shells included in bf function
  int NumKVecs; // number of k vectors included in bf function

  Vector<ComplexType> Rhok;

  Matrix<int> PairID;
  bool first;
  ///set of variables to be optimized
  opt_variables_type myVars;

  Backflow_ee_kSpace(ParticleSet& ions, ParticleSet& els): BackflowFunctionBase(ions,els)
  {
    first=true;
    Optimize=false;
    numParams=0;
    resize(NumTargets);
    NumGroups=els.groups();
    PairID.resize(NumTargets,NumTargets);
    for(int i=0; i<NumTargets; ++i)
      for(int j=0; j<NumTargets; ++j)
        PairID(i,j) = els.GroupID[i]*NumGroups+els.GroupID[j];
    offsetPrms.resize(NumGroups*NumGroups,0);
  }

  void initialize(ParticleSet&P, std::vector<RealType>& yk)
  {
    NumKShells = yk.size();
    Fk = yk;
    NumKVecs = P.SK->KLists.kshell[NumKShells+1];
    Rhok.resize(NumKVecs);
    if(Optimize)
      numParams = NumKShells;
  }

  void resize(int NT)
  {
    NumTargets=NT;
  }

  ~Backflow_ee_kSpace() {};

  void resetTargetParticleSet(ParticleSet& P)
  {
  }

  BackflowFunctionBase* makeClone(ParticleSet& tqp)
  {
    Backflow_ee_kSpace* clone = new Backflow_ee_kSpace(CenterSys,tqp);
    first=true;
    clone->resize(NumTargets);
//       clone->uniqueRadFun.resize(uniqueRadFun.size());
//       clone->RadFun.resize(RadFun.size());
    /*
    for(int i=0; i<uniqueRadFun.size(); i++)
      clone->uniqueRadFun[i] = new FT(*(uniqueRadFun[i]));
    for(int i=0; i<RadFun.size(); i++)
    {
      bool done=false;
      for(int k=0; k<uniqueRadFun.size(); k++)
        if(RadFun[i] == uniqueRadFun[k]) {
          done=true;
          clone->RadFun[i] = clone->uniqueRadFun[k];
          break;
        }
      if(!done) {
        APP_ABORT("Error cloning Backflow_ee_kSpace object. \n");
      }
    }
    */
    return clone;
  }

  void addFunc(int ia, int ib)
  {
    /*
          uniqueRadFun.push_back(rf);
          if(first) {
            // initialize all with rf the first time
            for(int i=0; i<RadFun.size(); i++)
              RadFun[i]=rf;
            first=false;
          } else {
            RadFun[ia*NumGroups+ib] = rf;
            RadFun[ib*NumGroups+ia] = rf;
          }
    */
  }

  void registerData(WFBufferType& buf)
  {
    /*
          FirstOfU = &(UIJ(0,0)[0]);
          LastOfU = FirstOfU + OHMMS_DIM*NumTargets*NumTargets;
          FirstOfA = &(AIJ(0,0)[0]);
          LastOfA = FirstOfA + OHMMS_DIM*OHMMS_DIM*NumTargets*NumTargets;
          FirstOfB = &(BIJ(0,0)[0]);
          LastOfB = FirstOfB + OHMMS_DIM*NumTargets*NumTargets;
          buf.add(FirstOfU,LastOfU);
          buf.add(FirstOfA,LastOfA);
          buf.add(FirstOfB,LastOfB);
    */
  }

  void reportStatus(std::ostream& os)
  {
    myVars.print(os);
  }

  void resetParameters(const opt_variables_type& active)
  {
    if(Optimize)
    {
      for (int i=0; i<Fk.size(); ++i)
      {
        int loc=myVars.where(i);
        if (loc>=0)
          Fk[i]=myVars[i]=active[loc];
      }
    }
  }

  void checkInVariables(opt_variables_type& active)
  {
    if(Optimize)
      active.insertFrom(myVars);
  }

  void checkOutVariables(const opt_variables_type& active)
  {
    if(Optimize)
      myVars.getIndex(active);
  }

  inline bool isOptimizable()
  {
    return Optimize;
  }

  inline int
  indexOffset()
  {
    if(Optimize)
      return myVars.where(0);
    else
      return 0;
  }

  inline void
  acceptMove(int iat, int UpdateMode)
  {
    int num;
    switch(UpdateMode)
    {
    case ORB_PBYP_RATIO:
      num = UIJ.rows();
      for(int i=0; i<num; i++)
      {
//            UIJ(iat,i) = UIJ_temp(i);
//            UIJ(i,iat) = -1.0*UIJ_temp(i);
      }
      break;
    case ORB_PBYP_PARTIAL:
      num = UIJ.rows();
      for(int i=0; i<num; i++)
      {
//            UIJ(iat,i) = UIJ_temp(i);
//            UIJ(i,iat) = -1.0*UIJ_temp(i);
      }
      num = AIJ.rows();
      for(int i=0; i<num; i++)
      {
//            AIJ(iat,i) = AIJ_temp(i);
//            AIJ(i,iat) = AIJ_temp(i);
      }
      break;
    case ORB_PBYP_ALL:
      num = UIJ.rows();
      for(int i=0; i<num; i++)
      {
//            UIJ(iat,i) = UIJ_temp(i);
//            UIJ(i,iat) = -1.0*UIJ_temp(i);
      }
      num = AIJ.rows();
      for(int i=0; i<num; i++)
      {
//            AIJ(iat,i) = AIJ_temp(i);
//            AIJ(i,iat) = AIJ_temp(i);
      }
      num = BIJ.rows();
      for(int i=0; i<num; i++)
      {
//            BIJ(iat,i) = BIJ_temp(i);
//            BIJ(i,iat) = -1.0*BIJ_temp(i);
      }
      break;
    default:
      num = UIJ.rows();
      for(int i=0; i<num; i++)
      {
//            UIJ(iat,i) = UIJ_temp(i);
//            UIJ(i,iat) = -1.0*UIJ_temp(i);
      }
      num = AIJ.rows();
      for(int i=0; i<num; i++)
      {
//            AIJ(iat,i) = AIJ_temp(i);
//            AIJ(i,iat) = AIJ_temp(i);
      }
      num = BIJ.rows();
      for(int i=0; i<num; i++)
      {
//            BIJ(iat,i) = BIJ_temp(i);
//            BIJ(i,iat) = -1.0*BIJ_temp(i);
      }
      break;
    }
//      UIJ_temp=0.0;
//      AIJ_temp=0.0;
//      BIJ_temp=0.0;
  }

  inline void
  restore(int iat, int UpdateType)
  {
//      UIJ_temp=0.0;
//      AIJ_temp=0.0;
//      BIJ_temp=0.0;
  }

  /** calculate quasi-particle coordinates only
   */
  inline void
  evaluate(const ParticleSet& P, ParticleSet& QP)
  {
#if defined(USE_REAL_STRUCT_FACTOR)
    APP_ABORT("Backflow_ee_kSpace::evaluate");
#else
    //memcopy if necessary but this is not so critcal
    copy(P.SK->rhok[0],P.SK->rhok[0]+NumKVecs,Rhok.data());
    for(int spec1=1; spec1<NumGroups; spec1++)
      accumulate_elements(P.SK->rhok[spec1],P.SK->rhok[spec1]+NumKVecs,Rhok.data());
    const KContainer::VContainer_t& Kcart(P.SK->KLists.kpts_cart);
    std::vector<int>& kshell(P.SK->KLists.kshell);
    for(int iel=0; iel<NumTargets; iel++)
    {
      const ComplexType* restrict eikr_ptr(P.SK->eikr[iel]);
      const ComplexType* restrict rhok_ptr(Rhok.data());
      int ki=0;
      for(int ks=0; ks<NumKShells; ks++)
      {
// don't understand factor of 2, ask Markus!!!
        RealType pre = 2.0*Fk[ks];
        for(; ki<kshell[ks+1]; ki++,eikr_ptr++,rhok_ptr++)
        {
          RealType ii=((*eikr_ptr).real()*(*rhok_ptr).imag()-(*eikr_ptr).imag()*(*rhok_ptr).real());
          QP.R[iel] -= pre*Kcart[ki]*ii;
        }
      }
    }
#endif
  }

  inline void
  evaluate(const ParticleSet& P, ParticleSet& QP, GradVector_t& Bmat, HessMatrix_t& Amat)
  {
    APP_ABORT("This shouldn't be called: Backflow_ee_kSpace::evaluate(Bmat)");
  }


  /** calculate quasi-particle coordinates, Bmat and Amat
   */
  inline void
  evaluate(const ParticleSet& P, ParticleSet& QP, GradMatrix_t& Bmat_full, HessMatrix_t& Amat)
  {
#if defined(USE_REAL_STRUCT_FACTOR)
    APP_ABORT("Backflow_ee_kSpace::evaluate");
#else
    //memcopy if necessary but this is not so critcal
    copy(P.SK->rhok[0],P.SK->rhok[0]+NumKVecs,Rhok.data());
    for(int spec1=1; spec1<NumGroups; spec1++)
      accumulate_elements(P.SK->rhok[spec1],P.SK->rhok[spec1]+NumKVecs,Rhok.data());
    const KContainer::VContainer_t& Kcart(P.SK->KLists.kpts_cart);
    std::vector<int>& kshell(P.SK->KLists.kshell);
    GradType fact;
    HessType kakb;
    for(int iel=0; iel<NumTargets; iel++)
    {
      const ComplexType* restrict eikr_ptr(P.SK->eikr[iel]);
      const ComplexType* restrict rhok_ptr(Rhok.data());
      const RealType* restrict ksq_ptr(P.SK->KLists.ksq.data());
      int ki=0;
      for(int ks=0; ks<NumKShells; ks++,ksq_ptr++)
      {
// don't understand factor of 2, ask Markus!!!
        RealType pre = 2.0*Fk[ks];
        RealType k2 = *ksq_ptr;
        for(; ki<kshell[ks+1]; ki++,eikr_ptr++,rhok_ptr++)
        {
          RealType rr=((*eikr_ptr).real()*(*rhok_ptr).real()+(*eikr_ptr).imag()*(*rhok_ptr).imag());
          RealType ii=((*eikr_ptr).real()*(*rhok_ptr).imag()-(*eikr_ptr).imag()*(*rhok_ptr).real());
          // quasiparticle
          QP.R[iel] -= pre*Kcart[ki]*ii;
          // B matrix
          convert(k2*pre*Kcart[ki],fact);
          convert(pre*outerProduct(Kcart[ki],Kcart[ki]),kakb);
          Bmat_full(iel,iel) += fact*(ii-1.0);
          // I can use symmetry to do only i<j, and then symmetrize
          // A matrix
          Amat(iel,iel) += kakb*(rr-1.0);
          // I can use symmetry to do only i<j, and then symmetrize
          for(int j=0; j<iel; j++)
          {
            ComplexType& eikrj = P.SK->eikr[j][ki];
            Bmat_full(iel,j) += fact*((*eikr_ptr).real()*(eikrj).imag()-(*eikr_ptr).imag()*(eikrj).real());
            Amat(iel,j) -= kakb*((*eikr_ptr).real()*(eikrj).real()+(*eikr_ptr).imag()*(eikrj).imag());
          }
          for(int j=iel+1; j<NumTargets; j++)
          {
            ComplexType& eikrj = P.SK->eikr[j][ki];
            Bmat_full(iel,j) += fact*((*eikr_ptr).real()*(eikrj).imag()-(*eikr_ptr).imag()*(eikrj).real());
            Amat(iel,j) -= kakb*((*eikr_ptr).real()*(eikrj).real()+(*eikr_ptr).imag()*(eikrj).imag());
          }
        }
      }
    }
#endif
  }

  /** calculate quasi-particle coordinates after pbyp move
   */
  inline void
  evaluatePbyP(const ParticleSet& P, ParticleSet::ParticlePos_t& newQP
               ,const std::vector<int>& index)
  {
  }

  /** calculate quasi-particle coordinates after pbyp move
   */
  inline void
  evaluatePbyP(const ParticleSet& P, int iat, ParticleSet::ParticlePos_t& newQP)
  {
  }

  /** calculate quasi-particle coordinates and Amat after pbyp move
   */
  inline void
  evaluatePbyP(const ParticleSet& P, ParticleSet::ParticlePos_t& newQP
               ,const std::vector<int>& index, HessMatrix_t& Amat)
  {
  }

  /** calculate quasi-particle coordinates and Amat after pbyp move
   */
  inline void
  evaluatePbyP(const ParticleSet& P, int iat, ParticleSet::ParticlePos_t& newQP
               , HessMatrix_t& Amat)
  {
  }

  /** calculate quasi-particle coordinates and Amat after pbyp move
   */
  inline void
  evaluatePbyP(const ParticleSet& P, ParticleSet::ParticlePos_t& newQP
               ,const std::vector<int>& index, GradMatrix_t& Bmat, HessMatrix_t& Amat)
  {
  }

  /** calculate quasi-particle coordinates and Amat after pbyp move
   */
  inline void
  evaluatePbyP(const ParticleSet& P, int iat, ParticleSet::ParticlePos_t& newQP
               , GradMatrix_t& Bmat, HessMatrix_t& Amat)
  {
  }

  /** calculate only Bmat
   *  This is used in pbyp moves, in updateBuffer()
   */
  inline void
  evaluateBmatOnly(const ParticleSet& P,GradMatrix_t& Bmat_full)
  {
  }

  /** calculate quasi-particle coordinates, Bmat and Amat
   *  calculate derivatives wrt to variational parameters
   */
  inline void
  evaluateWithDerivatives(const ParticleSet& P, ParticleSet& QP, GradMatrix_t& Bmat_full, HessMatrix_t& Amat, GradMatrix_t& Cmat, GradMatrix_t& Ymat, HessArray_t& Xmat)
  {
  }

};

}

#endif
