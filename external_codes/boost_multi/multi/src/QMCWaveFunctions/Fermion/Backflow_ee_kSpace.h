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
    NumKVecs   = P.getSimulationCell().getKLists().kshell[NumKShells + 1];
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
    throw std::runtime_error("Backflow_ee_kSpace::evaluate not implemented. There was an implementation with"
                             " complex-valued storage that may be resurrected using real-valued storage.");
  }

  inline void evaluate(const ParticleSet& P, ParticleSet& QP, GradVector& Bmat, HessMatrix& Amat)
  {
    APP_ABORT("This shouldn't be called: Backflow_ee_kSpace::evaluate(Bmat)");
  }


  /** calculate quasi-particle coordinates, Bmat and Amat
   */
  inline void evaluate(const ParticleSet& P, ParticleSet& QP, GradMatrix& Bmat_full, HessMatrix& Amat) override
  {
    throw std::runtime_error("Backflow_ee_kSpace::evaluate not implemented. There was an implementation with"
                             " complex-valued storage that may be resurrected using real-valued storage.");
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
