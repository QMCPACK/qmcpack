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
    
    
#ifndef QMCPLUSPLUS_BACKFLOW_TRANSFORMATION_H
#define QMCPLUSPLUS_BACKFLOW_TRANSFORMATION_H
#include "Particle/MCWalkerConfiguration.h"
#include "Utilities/ProgressReportEngine.h"
#include "OhmmsData/AttributeSet.h"
#include "QMCWaveFunctions/TrialWaveFunction.h"
#include "QMCWaveFunctions/OrbitalSetTraits.h"
#include "QMCWaveFunctions/Fermion/BackflowFunctionBase.h"
#include "QMCWaveFunctions/Fermion/Backflow_ee.h"
#include "QMCWaveFunctions/Fermion/Backflow_eI.h"
#include "QMCWaveFunctions/Fermion/GaussianFunctor.h"
#include "QMCWaveFunctions/Jastrow/BsplineFunctor.h"
#include "Particle/ParticleSet.h"
#include "Configuration.h"
#include <map>
#include <cmath>
#include "OhmmsPETE/OhmmsArray.h"

namespace qmcplusplus
{

class BackflowTransformation  //: public OrbitalSetTraits<QMCTraits::ValueType>
{

public:

  typedef BackflowFunctionBase::WFBufferType WFBufferType;

  // All BF quantities should be real, so eliminating complex (ValueType) possibility
  enum {DIM=OHMMS_DIM};
  typedef OHMMS_PRECISION                RealType;
  typedef int                            IndexType;
  typedef TinyVector<RealType,DIM>       PosType;
  typedef TinyVector<RealType,DIM>       GradType;
  typedef Tensor<RealType,DIM>           HessType;
  typedef Vector<IndexType>     IndexVector_t;
  typedef Vector<GradType>      GradVector_t;
  typedef Matrix<GradType>      GradMatrix_t;
  typedef Vector<HessType>      HessVector_t;
  typedef Matrix<HessType>      HessMatrix_t;

  typedef Array<HessType,3>       HessArray_t;

  typedef MCWalkerConfiguration::Walker_t Walker_t;
  typedef std::map<std::string,ParticleSet*>   PtclPoolType;
  //typedef Array<GradType,3>       GradArray_t;
  //typedef Array<PosType,3>        PosArray_t;

  ///number of quantum particles
  int NumTargets;

  /// active particle in pbyp moves
  int activeParticle;

  /// quasiparticle coordinates
  ParticleSet QP;

  /// Distance Table
  DistanceTableData* myTable;

  // number of variational parameters
  int numParams;

  /** current update mode */
  int UpdateMode;

  /** enum for a update mode */
  enum
  {
    ORB_PBYP_RATIO,   /*!< particle-by-particle ratio only */
    ORB_PBYP_ALL,     /*!< particle-by-particle, update Value-Gradient-Laplacian */
    ORB_PBYP_PARTIAL, /*!< particle-by-particle, update Value and Grdient */
    ORB_WALKER,    /*!< walker update */
    ORB_ALLWALKER  /*!< all walkers update */
  };

  // map index of variables from local arrays to outside world
  std::map<int,int> optIndexMap;

  // cutoff of radial funtions
  RealType cutOff;

  // pos of first optimizable variable in global array
  int numVarBefore;

  ParticleSet& targetPtcl;
//    PtclPoolType& ptclPool;

  // matrix of laplacians
  // /vec{B(i)} = sum_{k} /grad_{k}^2 /vec{x_i}
  GradVector_t Bmat;

  GradMatrix_t Bmat_full, Bmat_temp;

  // matrix of first derivatives
  // A(i,j)[a,b] = (Grad_i)_a (x_j)_b
  //               i,j:particle index
  //               a,b=(x,y,z)
// notice that A(i,j) is a symmetric matrix, improve later
  HessMatrix_t Amat, Amat_temp;

  // \nabla_a A_{i,j}^{\alpha,\beta}
  // derivative of A matrix with respect to var. prms.
  HessArray_t Xmat;

  // \sum_i \nabla_a B_{i,j}^{\alpha}
  GradMatrix_t Ymat;

  // \nabla_a x_i^{\alpha}
  GradMatrix_t Cmat;

  RealType *FirstOfP, *LastOfP;
  RealType *FirstOfA, *LastOfA;
  RealType *FirstOfB, *LastOfB;
  RealType *FirstOfA_temp, *LastOfA_temp;
  RealType *FirstOfB_temp, *LastOfB_temp;

  // Identity
  HessType HESS_ID;
  HessType DummyHess;

  std::vector<BackflowFunctionBase*> bfFuns;

  std::map<std::string,int> sources;
  std::vector<std::string> names;

  /// new qp coordinates for pbyp moves.
  ParticleSet::ParticlePos_t newQP;
  ParticleSet::ParticlePos_t oldQP;

  //Vector<PosType> storeQP;
  Vector<PosType> storeQP;

  /// store index of qp coordinates that changed during pbyp move
  std::vector<int> indexQP, index;

  opt_variables_type myVars;

//    BackflowTransformation(ParticleSet& els, PtclPoolType& pool):
//      targetPtcl(els),QP(els),/*ptclPool(pool),*/cutOff(0.0) {
//      myTable = DistanceTable::add(els,els);
//      NumTargets=els.getTotalNum();
//      Bmat.resize(NumTargets);
//      Bmat_full.resize(NumTargets,NumTargets);
//      Amat.resize(NumTargets,NumTargets);
//      newQP.resize(NumTargets);
//      oldQP.resize(NumTargets);
//      indexQP.resize(NumTargets);
//      HESS_ID.diagonal(1.0);
//      DummyHess=0.0;
//      numVarBefore=0;
//    }
  BackflowTransformation(ParticleSet& els):
    targetPtcl(els),QP(els),cutOff(0.0)
  {
    myTable = DistanceTable::add(els,DT_AOS);
    NumTargets=els.getTotalNum();
    Bmat.resize(NumTargets);
    Bmat_full.resize(NumTargets,NumTargets);
    Amat.resize(NumTargets,NumTargets);
    newQP.resize(NumTargets);
    oldQP.resize(NumTargets);
    indexQP.resize(NumTargets);
    HESS_ID.diagonal(1.0);
    DummyHess=0.0;
    numVarBefore=0;
  }

  void copyFrom(BackflowTransformation &tr)
  {
    cutOff=tr.cutOff;
    numParams=tr.numParams;
    numVarBefore=tr.numVarBefore;
    optIndexMap=tr.optIndexMap;
    bfFuns.resize((tr.bfFuns).size());
    std::vector<BackflowFunctionBase*>::iterator it((tr.bfFuns).begin());
    for(int i=0; i<(tr.bfFuns).size() ; i++,it++)
      bfFuns[i] = (*it)->makeClone(targetPtcl);
  }

// FIX FIX FIX
  BackflowTransformation* makeClone(ParticleSet& tqp)
  {
    BackflowTransformation *clone = new BackflowTransformation(tqp);
    clone->copyFrom(*this);
//       std::vector<BackflowFunctionBase*>::iterator it((bfFuns).begin());
//       for(int i=0; i<(bfFuns).size() ; i++,it++)
//       {
//         clone->bfFuns[i]->reportStatus(cerr);
//       }
    return clone;
  }

  ~BackflowTransformation() {};

  bool put(xmlNodePtr cur)
  {
    return true;
  }

  inline void
  acceptMove(const ParticleSet& P, int iat)
  {
    // update QP table
    // may be faster if I do this one qp at a time, for now do full update
    for(int i=0; i<NumTargets; i++)
      QP.R[i] = newQP[i];
    QP.update(0);
    indexQP.clear();
    switch(UpdateMode)
    {
    case ORB_PBYP_RATIO:
      break;
    case ORB_PBYP_PARTIAL:
      std::copy(FirstOfA_temp,LastOfA_temp,FirstOfA);
      break;
    case ORB_PBYP_ALL:
      std::copy(FirstOfA_temp,LastOfA_temp,FirstOfA);
      std::copy(FirstOfB_temp,LastOfB_temp,FirstOfB);
      break;
    default:
      std::copy(FirstOfA_temp,LastOfA_temp,FirstOfA);
      std::copy(FirstOfB_temp,LastOfB_temp,FirstOfB);
      break;
    }
    for(int i=0; i<bfFuns.size(); i++)
      bfFuns[i]->acceptMove(iat,UpdateMode);
  }

  inline void
  restore(int iat=0)
  {
    indexQP.clear();
    for(int i=0; i<bfFuns.size(); i++)
      bfFuns[i]->restore(iat,UpdateMode);
  }

  inline void checkInVariables(opt_variables_type& active)
  {
    for(int i=0; i<bfFuns.size(); i++)
      bfFuns[i]->checkInVariables(active);
  }

  inline void reportStatus(std::ostream& os)
  {
    for(int i=0; i<bfFuns.size(); i++)
      bfFuns[i]->reportStatus(os);
  }

  inline void checkOutVariables(const opt_variables_type& active)
  {
    for(int i=0; i<bfFuns.size(); i++)
      bfFuns[i]->checkOutVariables(active);
  }

  inline bool isOptimizable()
  {
    for(int i=0; i<bfFuns.size(); i++)
      if(bfFuns[i]->isOptimizable())
        return true;
    return false;
  }

  /** reset the distance table with a new target P
   */
  void resetTargetParticleSet(ParticleSet& P)
  {
    targetPtcl=P;
    myTable = DistanceTable::add(P,DT_AOS);
    for(int i=0; i<bfFuns.size(); i++)
      bfFuns[i]->resetTargetParticleSet(P);
  }


  void resetParameters(const opt_variables_type& active)
  {
    //reset each unique basis functions
    for(int i=0; i<bfFuns.size(); i++)
      if(bfFuns[i]->isOptimizable())
        bfFuns[i]->resetParameters(active);
  }

  void registerData(ParticleSet& P, WFBufferType& buf)
  {
    if(storeQP.size() == 0)
    {
      Bmat_temp.resize(NumTargets,NumTargets);
      Amat_temp.resize(NumTargets,NumTargets);
      storeQP.resize(NumTargets);
    }
    evaluate(P);
    FirstOfP = &(storeQP[0][0]);
    LastOfP = FirstOfP + OHMMS_DIM*NumTargets;
    FirstOfA = &(Amat(0,0)[0]);
    LastOfA = FirstOfA + OHMMS_DIM*OHMMS_DIM*NumTargets*NumTargets;
    FirstOfB = &(Bmat_full(0,0)[0]);
    LastOfB = FirstOfB + OHMMS_DIM*NumTargets*NumTargets;
    FirstOfA_temp = &(Amat_temp(0,0)[0]);
    LastOfA_temp = FirstOfA_temp + OHMMS_DIM*OHMMS_DIM*NumTargets*NumTargets;
    FirstOfB_temp = &(Bmat_temp(0,0)[0]);
    LastOfB_temp = FirstOfB_temp + OHMMS_DIM*NumTargets*NumTargets;
    for(int i=0; i<NumTargets; i++)
      storeQP[i] = QP.R[i];
    buf.add(FirstOfP,LastOfP);
    buf.add(FirstOfA,LastOfA);
    buf.add(FirstOfB,LastOfB);
    for(int i=0; i<bfFuns.size(); i++)
      bfFuns[i]->registerData(buf);
  }

  void updateBuffer(ParticleSet& P, WFBufferType& buf, bool redo)
  {
    //if(redo) evaluate(P);
    evaluate(P);
    for(int i=0; i<NumTargets; i++)
      storeQP[i] = QP.R[i];
    buf.put(FirstOfP,LastOfP);
    buf.put(FirstOfA,LastOfA);
    buf.put(FirstOfB,LastOfB);
    for(int i=0; i<bfFuns.size(); i++)
      bfFuns[i]->updateBuffer(buf);
  }

  void copyFromBuffer(ParticleSet& P, WFBufferType& buf)
  {
    buf.get(FirstOfP,LastOfP);
    buf.get(FirstOfA,LastOfA);
    buf.get(FirstOfB,LastOfB);
    for(int i=0; i<NumTargets; i++)
      QP.R[i] = storeQP[i];
    QP.update(0);
    for(int i=0; i<bfFuns.size(); i++)
      bfFuns[i]->copyFromBuffer(buf);
  }

  /** calculate quasi-particle coordinates only
   */
  inline void
  transformOnly(const ParticleSet& P)
  {
    for(int i=0; i<NumTargets; i++)
      QP.R[i] = P.R[i];
    for(int i=0; i<bfFuns.size(); i++)
      bfFuns[i]->evaluate(P,QP);
    QP.update(0);  // update distance tables
  }

  /** calculate new quasi-particle coordinates after pbyp move
   */
  inline void
  evaluatePbyP(const ParticleSet& P, int iat)
  //evaluatePbyP( ParticleSet& P, int iat)
  {
    UpdateMode = ORB_PBYP_RATIO;
    // there should be no need for this, but there is (missing calls in QMCHam...)
    for(int i=0; i<bfFuns.size(); i++)
      bfFuns[i]->restore(iat,UpdateMode);
    activeParticle=iat;
    for(int i=0; i<NumTargets; i++)
      oldQP[i] = newQP[i] = QP.R[i];
    newQP[iat] += myTable->Temp[iat].dr1;
    indexQP.clear();
    for(int i=0; i<bfFuns.size(); i++)
      bfFuns[i]->evaluatePbyP(P,iat,newQP);
    for(int jat=0; jat<NumTargets; jat++)
    {
      // make direct routine in OhmmsPETE later
      RealType dr = std::sqrt( dot( newQP[jat]-QP.R[jat], newQP[jat]-QP.R[jat] ) );
      if( dr > 1e-10 )
        indexQP.push_back(jat);
    }
//debug
    /*
    dummyQP2.R = P.R;
    dummyQP2.update();
    resetTargetParticleSet(dummyQP2);
    evaluate(P,dummyQP);
    resetTargetParticleSet(P);
    std::cout <<"index: ";
    for(int i=0; i<indexQP.size(); i++) std::cout <<indexQP[i] <<" ";
    std::cout << std::endl;
    for(int jat=0; jat<NumTargets; jat++)
      std::cout <<jat <<"  "
      <<(newQP[jat]-dummyQP.R[jat]) <<" " <<newQP[jat] <<" " <<QP.R[jat] <<"\n";
    for(int i=0; i<NumTargets; i++) newQP[i] = dummyQP.R[i];
    * /
    indexQP.clear();
    indexQP.push_back(iat); // set in the beggining by default
    for(int jat=0; jat<NumTargets; jat++) {
      if(jat!=iat) // && myTable->Temp[jat].r1 < cutOff )
        indexQP.push_back(jat);
    }
    */
  }

  /** calculate new quasi-particle coordinates after pbyp move
   */
  inline void
  evaluatePbyPWithGrad(const ParticleSet& P, int iat)
  {
    UpdateMode = ORB_PBYP_PARTIAL;
    // there should be no need for this, but there is (missing calls in QMCHam...)
    for(int i=0; i<bfFuns.size(); i++)
      bfFuns[i]->restore(iat,UpdateMode);
    activeParticle=iat;
    for(int i=0; i<NumTargets; i++)
      oldQP[i] = newQP[i] = QP.R[i];
    newQP[iat] += myTable->Temp[iat].dr1;
    indexQP.clear();
    std::copy(FirstOfA,LastOfA,FirstOfA_temp);
    for(int i=0; i<bfFuns.size(); i++)
      bfFuns[i]->evaluatePbyP(P,iat,newQP,Amat_temp);
    for(int jat=0; jat<NumTargets; jat++)
    {
      // make direct routine in OhmmsPETE later
      RealType dr = std::sqrt( dot( newQP[jat]-QP.R[jat], newQP[jat]-QP.R[jat] ) );
      if( dr > 1e-10 )
        indexQP.push_back(jat);
    }
  }

  /** calculate new quasi-particle coordinates after pbyp move
   */
  inline void
  evaluatePbyPAll(const ParticleSet& P, int iat)
  {
    UpdateMode = ORB_PBYP_ALL;
    // there should be no need for this, but there is (missing calls in QMCHam...)
    for(int i=0; i<bfFuns.size(); i++)
      bfFuns[i]->restore(iat,UpdateMode);
    activeParticle=iat;
    for(int i=0; i<NumTargets; i++)
      oldQP[i] = newQP[i] = QP.R[i];
    newQP[iat] += myTable->Temp[iat].dr1;
    indexQP.clear();
    std::copy(FirstOfA,LastOfA,FirstOfA_temp);
    std::copy(FirstOfB,LastOfB,FirstOfB_temp);
    for(int i=0; i<bfFuns.size(); i++)
      bfFuns[i]->evaluatePbyP(P,iat,newQP,Bmat_temp,Amat_temp);
    for(int jat=0; jat<NumTargets; jat++)
    {
      // make direct routine in OhmmsPETE later
      RealType dr = std::sqrt( dot( newQP[jat]-QP.R[jat], newQP[jat]-QP.R[jat] ) );
      if( dr > 1e-10 )
        indexQP.push_back(jat);
    }
  }


  /** calculate only Bmat. Assume that QP and Amat are current
   *  This is used in pbyp moves, in updateBuffer()
   */
  inline void
  evaluateBmatOnly(const ParticleSet& P, int iat)
  {
    Bmat_full=0.0;
    for(int i=0; i<bfFuns.size(); i++)
      bfFuns[i]->evaluateBmatOnly(P,Bmat_full);
  }

  /** calculate quasi-particle coordinates, Bmat and Amat
   */
  inline void
  evaluate(const ParticleSet& P)
  {
    Bmat=0.0;
    Amat=0.0;
    Bmat_full=0.0;
    QP.R=P.R;
    for(int i=0; i<NumTargets; i++)
    {
      //QP.R[i] = P.R[i];
      Amat(i,i).diagonal(1.0);
    }
    for(int i=0; i<bfFuns.size(); i++)
      bfFuns[i]->evaluate(P,QP,Bmat_full,Amat);
//      std::cerr <<"P.R \n";
//      std::cerr <<P.R[0] << std::endl;
//      std::cerr <<"QP.R " << std::endl;
//      std::cerr <<QP.R[0] << std::endl;
//      std::cerr <<omp_get_thread_num()<<" "<<P.R[0]-QP.R[0] << std::endl;
//      APP_ABORT("TESTING BF \n");
    /*Bmat=0.0;
    Amat=0.0;
    Bmat_full=0.0;
    for(int i=0; i<NumTargets; i++) {
      Amat(i,i).diagonal(1.0);
    }*/
    /*
          // testing bf
          for(int i=0; i<NumTargets; i++) {
            std::cout <<"i: " <<i << std::endl;
            std::cout <<P.R[i] << std::endl;
            std::cout <<QP.R[i] << std::endl;
            std::cout <<P.R[i]-QP.R[i] << std::endl;
          }
          //
    */
    QP.update(0);  // update distance tables
  }

  /** calculate quasi-particle coordinates and store in Pnew
   */
  inline void
  evaluate(const ParticleSet& P, ParticleSet& Pnew)
  {
    Pnew.R=P.R;
    for(int i=0; i<bfFuns.size(); i++)
      bfFuns[i]->evaluate(P,Pnew);
    Pnew.update(0);
  }

  inline void
  evaluateDerivatives(const ParticleSet& P)
  {
    if(Cmat.size() == 0)
      // initialize in the first call
    {
      // assumes that all BF parameters are packed together in
      // active variable set. is this always correct???
      numParams=0;
      for(int i=0; i<bfFuns.size(); i++)
      {
        int tmp = bfFuns[i]->setParamIndex(numParams);
        numParams+=tmp;
      }
      numVarBefore = bfFuns[0]->indexOffset();
      //app_log() <<"numVarBefore: " <<numVarBefore << std::endl;
      for(int i=0; i<numParams; i++)
      {
        optIndexMap[i] = i+numVarBefore;
        //app_log() <<"prm, map: " <<i <<"  " <<optIndexMap[i] << std::endl;
      }
      Cmat.resize(numParams,NumTargets);
      Xmat.resize(numParams,NumTargets,NumTargets);
      Ymat.resize(numParams,NumTargets);
    }
    // Uncomment to test calculation of Cmat,Xmat,Ymat
    //testDeriv(P);
    Bmat=0.0;
    Amat=0.0;
    Bmat_full=0.0;
    Cmat=0.0;
    Ymat=0.0;
    for(int i=0; i<Xmat.size(); i++)
      Xmat(i)=0;
    for(int i=0; i<NumTargets; i++)
    {
      QP.R[i] = P.R[i];
      Amat(i,i).diagonal(1.0);
    }
    for(int i=0; i<bfFuns.size(); i++)
      bfFuns[i]->evaluateWithDerivatives(P,QP,Bmat_full,Amat,Cmat,Ymat,Xmat);
    QP.update(0);
  }

  void testDeriv(const ParticleSet& P)
  {
    if(Cmat.size() == 0)
      // initialize in the first call
    {
      Cmat.resize(numParams,NumTargets);
      Xmat.resize(numParams,NumTargets,NumTargets);
      Ymat.resize(numParams,NumTargets);
    }
    Bmat=0.0;
    Amat=0.0;
    Bmat_full=0.0;
    Cmat=0.0;
    Ymat=0.0;
//       Xmat=DummyHess;
    for(int i=0; i<Xmat.size(); i++)
      Xmat(i)=0;
    for(int i=0; i<NumTargets; i++)
    {
      QP.R[i] = P.R[i];
      Amat(i,i).diagonal(1.0);
    }
    for(int i=0; i<bfFuns.size(); i++)
      bfFuns[i]->evaluateWithDerivatives(P,QP,Bmat_full,Amat,Cmat,Ymat,Xmat);
    ParticleSet::ParticlePos_t qp_0;
    ParticleSet::ParticlePos_t qp_1;
    ParticleSet::ParticlePos_t qp_2;
    GradMatrix_t Bmat_full_1;
    HessMatrix_t Amat_1;
    GradMatrix_t Bmat_full_2;
    HessMatrix_t Amat_2;
    RealType dh = 0.00001;
    qp_0.resize(NumTargets);
    qp_1.resize(NumTargets);
    qp_2.resize(NumTargets);
    Bmat_full_1.resize(NumTargets,NumTargets);
    Bmat_full_2.resize(NumTargets,NumTargets);
    Amat_1.resize(NumTargets,NumTargets);
    Amat_2.resize(NumTargets,NumTargets);
    for(int i=0; i<NumTargets; i++)
    {
      qp_0[i] = QP.R[i];
    }
    app_log() <<" Testing derivatives of backflow transformation. \n";
    app_log() <<" Numtargets: " <<NumTargets << std::endl;
    opt_variables_type wfVars,wfvar_prime;
    checkInVariables(wfVars);
    checkOutVariables(wfVars);
    int Nvars= wfVars.size();
    wfvar_prime= wfVars;
    wfVars.print(std::cout);
    for(int i=0; i<Nvars; i++)
    {
      for (int j=0; j<Nvars; j++)
        wfvar_prime[j]=wfVars[j];
      wfvar_prime[i] = wfVars[i]+ dh;
      resetParameters(wfvar_prime);
      Bmat_full_1=0.0;
      Amat_1=0.0;
      for(int k=0; k<NumTargets; k++)
      {
        QP.R[k] = P.R[k];
        Amat_1(k,k).diagonal(1.0);
      }
      for(int k=0; k<bfFuns.size(); k++)
        bfFuns[k]->evaluate(P,QP,Bmat_full_1,Amat_1);
      for(int k=0; k<NumTargets; k++)
        qp_1[k] = QP.R[k];
      for (int j=0; j<Nvars; j++)
        wfvar_prime[j]=wfVars[j];
      wfvar_prime[i] = wfVars[i] - dh;
      resetParameters(wfvar_prime);
      Bmat_full_2=0.0;
      Amat_2=0.0;
      for(int k=0; k<NumTargets; k++)
      {
        QP.R[k] = P.R[k];
        Amat_2(k,k).diagonal(1.0);
      }
      for(int k=0; k<bfFuns.size(); k++)
        bfFuns[k]->evaluate(P,QP,Bmat_full_2,Amat_2);
      for(int k=0; k<NumTargets; k++)
        qp_2[k] = QP.R[k];
      app_log() <<"Cmat: \n"
                <<"i, AvDiff, max: \n";
      //2011-07-17: what is the proper data type?
      RealType df,av=0.0,cnt=0.0;
      RealType maxD=-100.0;
      const RealType ConstOne(1.0);
      for(int k=0; k<NumTargets; k++)
      {
        for(int q=0; q<OHMMS_DIM; q++)
        {
          cnt += ConstOne;
          df=(( (qp_1[k])[q] - (qp_2[k])[q] )/(2.0*dh)-Cmat(i,k)[q]);
          av+=df;
          if( std::abs(df) > maxD )
            maxD=std::abs(df);
          //app_log() <<k <<"  " <<q <<"   "
          //          <<( (qp_1[k])[q] - (qp_2[k])[0] )/(2.0*dh)   <<"  "
          //          <<Cmat(i,k)[q] <<"  " <<(( (qp_1[k])[q] - (qp_2[k])[q] )/(2.0*dh)-Cmat(i,k)[q]) << std::endl;
        }
      }
      app_log() <<i <<"  " <<av/cnt <<"  " <<maxD << std::endl;
      av=cnt=maxD=0.0;
      app_log() <<"Ymat: \n";
      for(int k=0; k<NumTargets; k++)
      {
        for(int q=0; q<3; q++)
        {
          RealType dB=0.0;
          for(int j=0; j<NumTargets; j++)
            dB+= (Bmat_full_1(j,k)[q] - Bmat_full_2(j,k)[q]);
          cnt+=ConstOne;
          df=(dB/(2.0*dh)-Ymat(i,k)[q]);
          av+=df;
          if( std::abs(df) > maxD )
            maxD=std::abs(df);
          //app_log() <<k <<"  " <<q <<"   "
          //        <<dB/(2.0*dh)   <<"  "
          //        <<Ymat(i,k)[q] <<"  " <<(dB/(2.0*dh)-Ymat(i,k)[q]) << std::endl;
        }
      }
      app_log() <<i <<"  " <<av/cnt <<"  " <<maxD << std::endl;
      av=cnt=maxD=0.0;
      app_log() <<"Xmat: \n";
      for(int k1=0; k1<NumTargets; k1++)
        for(int k2=0; k2<NumTargets; k2++)
        {
          for(int q1=0; q1<3; q1++)
          {
            for(int q2=0; q2<3; q2++)
            {
              RealType dB=(Amat_1(k1,k2))(q1,q2) - (Amat_2(k1,k2))(q1,q2);
              cnt+=ConstOne;
              df=(dB/(2.0*dh)-(Xmat(i,k1,k2))(q1,q2));
              av+=df;
              if( std::abs(df) > maxD )
                maxD=std::abs(df);
              //app_log() <<k1 <<"  " <<k2 <<"  " <<q1 <<"  " <<q2 <<"   "
              //        <<(Xmat(i,k1,k2))(q1,q2) <<"  " <<(dB/(2.0*dh)-(Xmat(i,k1,k2))(q1,q2)) << std::endl;
            }
          }
        }
      app_log() <<i <<"  " <<av/cnt <<"  " <<maxD << std::endl;
      av=cnt=maxD=0.0;
    }
  }

  void testPbyP(ParticleSet& P)
  {
    GradMatrix_t Bmat_full_0;
    HessMatrix_t Amat_0;
    GradMatrix_t Bmat_full_1;
    HessMatrix_t Amat_1;
    ParticleSet::ParticlePos_t qp_0;
    ParticleSet::ParticlePos_t qp_1;
    ParticleSet::ParticlePos_t qp_2,qp_3;
    qp_0.resize(NumTargets);
    qp_1.resize(NumTargets);
    qp_2.resize(NumTargets);
    qp_3.resize(NumTargets);
    Bmat_full_0.resize(NumTargets,NumTargets);
    Bmat_full_1.resize(NumTargets,NumTargets);
    Amat_0.resize(NumTargets,NumTargets);
    Amat_1.resize(NumTargets,NumTargets);
    P.update();
    Walker_t::WFBuffer_t tbuffer;
    size_t BufferCursor=tbuffer.current();
    registerData(P,tbuffer);
    tbuffer.rewind(BufferCursor);
    updateBuffer(P,tbuffer,true);
    qp_3 = P.R;
    evaluate(P);
    qp_2 = QP.R;
    app_log() <<"after 1st eval: " <<cutOff << std::endl;
    for(int jat=0; jat<NumTargets; jat++)
      app_log() <<jat <<"  "
                <<P.R[jat]-QP.R[jat]  << std::endl;
    //for(int  iat=0; iat<NumTargets; iat++) {
    for(int  iat=0; iat<1; iat++)
    {
      PosType dr;
      dr[0]=0.1;
      dr[1]=0.05;
      dr[2]=-0.3;
      P.makeMove(iat,dr);
      app_log() <<"Move: " <<myTable->Temp[iat].dr1 << std::endl;
      app_log() <<"cutOff: " <<cutOff << std::endl;
      for(int jat=0; jat<NumTargets; jat++)
        app_log() <<jat <<"  " <<myTable->Temp[jat].r1 << std::endl;
      //evaluatePbyP(P,iat);
      evaluatePbyPWithGrad(P,iat);
      app_log() <<"Moving: ";
      for(int i=0; i<indexQP.size(); i++)
        app_log() <<indexQP[i] <<" ";
      app_log() << std::endl;
      acceptMove(P,iat);
      P.acceptMove(iat);
    }
    qp_0 = QP.R;
    Amat_0 = Amat;
    tbuffer.rewind(BufferCursor);
    updateBuffer(P,tbuffer,false);
    P.update();
    evaluate(P);
    Amat_1 = Amat_0 - Amat;
    qp_1 = QP.R - qp_0;
    RealType qpdiff = Dot(qp_1,qp_1);
    RealType Amdiff = 0.0;
    for(int i=0; i<NumTargets; i++)
      for(int k=0; k<NumTargets; k++)
        for(int j=0; j<OHMMS_DIM*OHMMS_DIM; j++)
          Amdiff += Amat_1(i,k)[j]*Amat_1(i,k)[j];
    app_log() <<"Error in pbyp QP transformation: " <<qpdiff << std::endl;
    app_log() <<"Error in pbyp QP Amat: " <<Amdiff << std::endl;
    app_log() <<"i, diff, newPbyP, newEval: \n";
    for(int i=0; i<NumTargets; i++)
      app_log() <<i <<"\n"
                <<qp_0[i]-QP.R[i] <<"\n"
                <<qp_0[i] <<"\n"
                <<QP.R[i] << std::endl << std::endl;
    APP_ABORT("Finished BackflowTransformation::testPbyP() \n.");
  }

};

}

#endif
