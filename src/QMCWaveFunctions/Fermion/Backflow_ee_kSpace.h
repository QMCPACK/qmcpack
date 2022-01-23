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
#include "LongRange/StructFact.h"
#include "Message/Communicate.h"
#include <cmath>

namespace qmcplusplus
{
class Backflow_ee_kSpace : public BackflowFunctionBase
{
  using ComplexType = QMCTraits::ComplexType;
  ///typedef for real values
  //using real_type = optimize::VariableSet::real_type;
  ///typedef for variableset: this is going to be replaced
  using opt_variables_type = optimize::VariableSet;
  ///typedef for name-value lists
  using variable_map_type = optimize::VariableSet::variable_map_type;

public:
  //number of groups of the target particleset
  bool Optimize;
  int numParams;
  std::vector<RealType> Fk;
  std::vector<int> offsetPrms;
  int NumGroups;
  int NumKShells; // number of k shells included in bf function
  int NumKVecs;   // number of k vectors included in bf function

  Vector<ComplexType> Rhok;

  Matrix<int> PairID;
  ///set of variables to be optimized
  opt_variables_type myVars;

  Backflow_ee_kSpace(ParticleSet& ions, ParticleSet& els) : BackflowFunctionBase(ions, els)
  {
    Optimize  = false;
    numParams = 0;
    resize(NumTargets);
    NumGroups = els.groups();
    PairID.resize(NumTargets, NumTargets);
    for (int i = 0; i < NumTargets; ++i)
      for (int j = 0; j < NumTargets; ++j)
        PairID(i, j) = els.GroupID[i] * NumGroups + els.GroupID[j];
    offsetPrms.resize(NumGroups * NumGroups, 0);
  }

  void initialize(ParticleSet& P, std::vector<RealType>& yk)
  {
    NumKShells = yk.size();
    Fk         = yk;
    NumKVecs   = P.getSK().getKLists().kshell[NumKShells + 1];
    Rhok.resize(NumKVecs);
    if (Optimize)
      numParams = NumKShells;
  }

  void resize(int NT) { NumTargets = NT; }

  std::unique_ptr<BackflowFunctionBase> makeClone(ParticleSet& tqp) const override
  {
    auto clone = std::make_unique<Backflow_ee_kSpace>(CenterSys, tqp);
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

  void registerData(WFBufferType& buf) override
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

  void reportStatus(std::ostream& os) override { myVars.print(os); }

  void resetParameters(const opt_variables_type& active) override
  {
    if (Optimize)
    {
      for (int i = 0; i < Fk.size(); ++i)
      {
        int loc = myVars.where(i);
        if (loc >= 0)
        {
          myVars[i] = active[loc];
          Fk[i]     = std::real(myVars[i]);
          //Fk[i] = myVars[i] = active[loc];
        }
      }
    }
  }

  void checkInVariables(opt_variables_type& active) override
  {
    if (Optimize)
      active.insertFrom(myVars);
  }

  void checkOutVariables(const opt_variables_type& active) override
  {
    if (Optimize)
      myVars.getIndex(active);
  }

  inline bool isOptimizable() override { return Optimize; }

  inline int indexOffset() override
  {
    if (Optimize)
      return myVars.where(0);
    else
      return 0;
  }

  inline void acceptMove(int iat, int UpdateMode) override
  {
    int num;
    switch (UpdateMode)
    {
    case ORB_PBYP_RATIO:
      num = UIJ.rows();
      for (int i = 0; i < num; i++)
      {
        //            UIJ(iat,i) = UIJ_temp(i);
        //            UIJ(i,iat) = -1.0*UIJ_temp(i);
      }
      break;
    case ORB_PBYP_PARTIAL:
      num = UIJ.rows();
      for (int i = 0; i < num; i++)
      {
        //            UIJ(iat,i) = UIJ_temp(i);
        //            UIJ(i,iat) = -1.0*UIJ_temp(i);
      }
      num = AIJ.rows();
      for (int i = 0; i < num; i++)
      {
        //            AIJ(iat,i) = AIJ_temp(i);
        //            AIJ(i,iat) = AIJ_temp(i);
      }
      break;
    case ORB_PBYP_ALL:
      num = UIJ.rows();
      for (int i = 0; i < num; i++)
      {
        //            UIJ(iat,i) = UIJ_temp(i);
        //            UIJ(i,iat) = -1.0*UIJ_temp(i);
      }
      num = AIJ.rows();
      for (int i = 0; i < num; i++)
      {
        //            AIJ(iat,i) = AIJ_temp(i);
        //            AIJ(i,iat) = AIJ_temp(i);
      }
      num = BIJ.rows();
      for (int i = 0; i < num; i++)
      {
        //            BIJ(iat,i) = BIJ_temp(i);
        //            BIJ(i,iat) = -1.0*BIJ_temp(i);
      }
      break;
    default:
      num = UIJ.rows();
      for (int i = 0; i < num; i++)
      {
        //            UIJ(iat,i) = UIJ_temp(i);
        //            UIJ(i,iat) = -1.0*UIJ_temp(i);
      }
      num = AIJ.rows();
      for (int i = 0; i < num; i++)
      {
        //            AIJ(iat,i) = AIJ_temp(i);
        //            AIJ(i,iat) = AIJ_temp(i);
      }
      num = BIJ.rows();
      for (int i = 0; i < num; i++)
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

  inline void restore(int iat, int UpdateType) override
  {
    //      UIJ_temp=0.0;
    //      AIJ_temp=0.0;
    //      BIJ_temp=0.0;
  }

  /** calculate quasi-particle coordinates only
   */
  inline void evaluate(const ParticleSet& P, ParticleSet& QP) override
  {
#if defined(USE_REAL_STRUCT_FACTOR)
    APP_ABORT("Backflow_ee_kSpace::evaluate");
#else
    //memcopy if necessary but this is not so critcal
    copy(P.getSK().rhok[0], P.getSK().rhok[0] + NumKVecs, Rhok.data());
    for (int spec1 = 1; spec1 < NumGroups; spec1++)
      accumulate_elements(P.getSK().rhok[spec1], P.getSK().rhok[spec1] + NumKVecs, Rhok.data());
    const auto& Kcart(P.getSK().getKLists().kpts_cart);
    std::vector<int>& kshell(P.getSK().getKLists().kshell);
    for (int iel = 0; iel < NumTargets; iel++)
    {
      const ComplexType* restrict eikr_ptr(P.getSK().eikr[iel]);
      const ComplexType* restrict rhok_ptr(Rhok.data());
      int ki = 0;
      for (int ks = 0; ks < NumKShells; ks++)
      {
        // don't understand factor of 2, ask Markus!!!
        RealType pre = 2.0 * Fk[ks];
        for (; ki < kshell[ks + 1]; ki++, eikr_ptr++, rhok_ptr++)
        {
          RealType ii = ((*eikr_ptr).real() * (*rhok_ptr).imag() - (*eikr_ptr).imag() * (*rhok_ptr).real());
          QP.R[iel] -= pre * Kcart[ki] * ii;
        }
      }
    }
#endif
  }

  inline void evaluate(const ParticleSet& P, ParticleSet& QP, GradVector& Bmat, HessMatrix& Amat)
  {
    APP_ABORT("This shouldn't be called: Backflow_ee_kSpace::evaluate(Bmat)");
  }


  /** calculate quasi-particle coordinates, Bmat and Amat
   */
  inline void evaluate(const ParticleSet& P, ParticleSet& QP, GradMatrix& Bmat_full, HessMatrix& Amat) override
  {
#if defined(USE_REAL_STRUCT_FACTOR)
    APP_ABORT("Backflow_ee_kSpace::evaluate");
#else
    //memcopy if necessary but this is not so critcal
    copy(P.getSK().rhok[0], P.getSK().rhok[0] + NumKVecs, Rhok.data());
    for (int spec1 = 1; spec1 < NumGroups; spec1++)
      accumulate_elements(P.getSK().rhok[spec1], P.getSK().rhok[spec1] + NumKVecs, Rhok.data());
    const auto& Kcart(P.getSK().getKLists().kpts_cart);
    std::vector<int>& kshell(P.getSK().getKLists().kshell);
    GradType fact;
    HessType kakb;
    for (int iel = 0; iel < NumTargets; iel++)
    {
      const ComplexType* restrict eikr_ptr(P.getSK().eikr[iel]);
      const ComplexType* restrict rhok_ptr(Rhok.data());
      const RealType* restrict ksq_ptr(P.getSK().getKLists().ksq.data());
      int ki = 0;
      for (int ks = 0; ks < NumKShells; ks++, ksq_ptr++)
      {
        // don't understand factor of 2, ask Markus!!!
        RealType pre = 2.0 * Fk[ks];
        RealType k2  = *ksq_ptr;
        for (; ki < kshell[ks + 1]; ki++, eikr_ptr++, rhok_ptr++)
        {
          RealType rr = ((*eikr_ptr).real() * (*rhok_ptr).real() + (*eikr_ptr).imag() * (*rhok_ptr).imag());
          RealType ii = ((*eikr_ptr).real() * (*rhok_ptr).imag() - (*eikr_ptr).imag() * (*rhok_ptr).real());
          // quasiparticle
          QP.R[iel] -= pre * Kcart[ki] * ii;
          // B matrix
          convert(k2 * pre * Kcart[ki], fact);
          convert(pre * outerProduct(Kcart[ki], Kcart[ki]), kakb);
          Bmat_full(iel, iel) += fact * (ii - 1.0);
          // I can use symmetry to do only i<j, and then symmetrize
          // A matrix
          Amat(iel, iel) += kakb * (rr - 1.0);
          // I can use symmetry to do only i<j, and then symmetrize
          for (int j = 0; j < iel; j++)
          {
            ComplexType& eikrj = P.getSK().eikr[j][ki];
            Bmat_full(iel, j) += fact * ((*eikr_ptr).real() * (eikrj).imag() - (*eikr_ptr).imag() * (eikrj).real());
            Amat(iel, j) -= kakb * ((*eikr_ptr).real() * (eikrj).real() + (*eikr_ptr).imag() * (eikrj).imag());
          }
          for (int j = iel + 1; j < NumTargets; j++)
          {
            ComplexType& eikrj = P.getSK().eikr[j][ki];
            Bmat_full(iel, j) += fact * ((*eikr_ptr).real() * (eikrj).imag() - (*eikr_ptr).imag() * (eikrj).real());
            Amat(iel, j) -= kakb * ((*eikr_ptr).real() * (eikrj).real() + (*eikr_ptr).imag() * (eikrj).imag());
          }
        }
      }
    }
#endif
  }

  /** calculate quasi-particle coordinates after pbyp move
   */
  inline void evaluatePbyP(const ParticleSet& P,
                           ParticleSet::ParticlePos& newQP,
                           const std::vector<int>& index) override
  {}

  /** calculate quasi-particle coordinates after pbyp move
   */
  inline void evaluatePbyP(const ParticleSet& P, int iat, ParticleSet::ParticlePos& newQP) override {}

  /** calculate quasi-particle coordinates and Amat after pbyp move
   */
  inline void evaluatePbyP(const ParticleSet& P,
                           ParticleSet::ParticlePos& newQP,
                           const std::vector<int>& index,
                           HessMatrix& Amat) override
  {}

  /** calculate quasi-particle coordinates and Amat after pbyp move
   */
  inline void evaluatePbyP(const ParticleSet& P, int iat, ParticleSet::ParticlePos& newQP, HessMatrix& Amat) override {}

  /** calculate quasi-particle coordinates and Amat after pbyp move
   */
  inline void evaluatePbyP(const ParticleSet& P,
                           ParticleSet::ParticlePos& newQP,
                           const std::vector<int>& index,
                           GradMatrix& Bmat,
                           HessMatrix& Amat) override
  {}

  /** calculate quasi-particle coordinates and Amat after pbyp move
   */
  inline void evaluatePbyP(const ParticleSet& P,
                           int iat,
                           ParticleSet::ParticlePos& newQP,
                           GradMatrix& Bmat,
                           HessMatrix& Amat) override
  {}

  /** calculate only Bmat
   *  This is used in pbyp moves, in updateBuffer()
   */
  inline void evaluateBmatOnly(const ParticleSet& P, GradMatrix& Bmat_full) override {}

  /** calculate quasi-particle coordinates, Bmat and Amat
   *  calculate derivatives wrt to variational parameters
   */
  inline void evaluateWithDerivatives(const ParticleSet& P,
                                      ParticleSet& QP,
                                      GradMatrix& Bmat_full,
                                      HessMatrix& Amat,
                                      GradMatrix& Cmat,
                                      GradMatrix& Ymat,
                                      HessArray& Xmat) override
  {}
};

} // namespace qmcplusplus

#endif
