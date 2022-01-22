//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Ken Esler, kpesler@gmail.com, University of Illinois at Urbana-Champaign
//                    Miguel Morales, moralessilva2@llnl.gov, Lawrence Livermore National Laboratory
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Ken Esler, kpesler@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////


#include "QMCLinearOptimize.h"
#include "Particle/HDFWalkerIO.h"
#include "OhmmsData/AttributeSet.h"
#include "Message/CommOperators.h"

#include "QMCDrivers/VMC/VMC.h"
#include "QMCDrivers/WFOpt/QMCCostFunction.h"

//#include "QMCDrivers/VMC/VMCSingle.h"
//#include "QMCDrivers/QMCCostFunctionSingle.h"
#include "QMCHamiltonians/HamiltonianPool.h"
#include "CPU/Blasf.h"
#include "Numerics/MatrixOperators.h"
#include <cassert>
#if defined(QMC_CUDA)
#include "QMCDrivers/VMC/VMC_CUDA.h"
#include "QMCDrivers/WFOpt/QMCCostFunctionCUDA.h"
#endif
#include "Numerics/LinearFit.h"
#include <iostream>
#include <fstream>

/*#include "Message/Communicate.h"*/

namespace qmcplusplus
{
QMCLinearOptimize::QMCLinearOptimize(MCWalkerConfiguration& w,
                                     TrialWaveFunction& psi,
                                     QMCHamiltonian& h,
                                     Communicate* comm,
                                     const std::string& QMC_driver_type)
    : QMCDriver(w, psi, h, comm, QMC_driver_type),
      PartID(0),
      NumParts(1),
      wfNode(NULL),
      optNode(NULL),
      param_tol(1e-4),
      generate_samples_timer_(*timer_manager.createTimer("QMCLinearOptimize::GenerateSamples", timer_level_medium)),
      initialize_timer_(*timer_manager.createTimer("QMCLinearOptimize::Initialize", timer_level_medium)),
      eigenvalue_timer_(*timer_manager.createTimer("QMCLinearOptimize::Eigenvalue", timer_level_medium)),
      line_min_timer_(*timer_manager.createTimer("QMCLinearOptimize::Line_Minimization", timer_level_medium)),
      cost_function_timer_(*timer_manager.createTimer("QMCLinearOptimize::CostFunction", timer_level_medium))
{
  IsQMCDriver = false;
  //     //set the optimization flag
  qmc_driver_mode.set(QMC_OPTIMIZE, 1);
  //read to use vmc output (just in case)
  m_param.add(param_tol, "alloweddifference");
  //Set parameters for line minimization:
}

/** Add configuration files for the optimization
* @param a root of a hdf5 configuration file
*/
void QMCLinearOptimize::addConfiguration(const std::string& a)
{
  if (a.size())
    ConfigFile.push_back(a);
}

void QMCLinearOptimize::start()
{
  {
    //generate samples
    ScopedTimer local(generate_samples_timer_);
    generateSamples();
    //store active number of walkers
    NumOfVMCWalkers = W.getActiveWalkers();
  }

  app_log() << "<opt stage=\"setup\">" << std::endl;
  app_log() << "  <log>" << std::endl;
  //reset the rootname
  optTarget->setRootName(RootName);
  optTarget->setWaveFunctionNode(wfNode);
  app_log() << "   Reading configurations from h5FileRoot " << h5FileRoot << std::endl;
  {
    //get configuration from the previous run
    ScopedTimer local(initialize_timer_);
    Timer t2;
    optTarget->getConfigurations(h5FileRoot);
    optTarget->setRng(vmcEngine->getRngRefs());
    optTarget->checkConfigurations();
    // check recomputed variance against VMC
    auto sigma2_vmc   = vmcEngine->getBranchEngine()->vParam[SimpleFixedNodeBranch::SBVP::SIGMA2];
    auto sigma2_check = optTarget->getVariance();
    if (sigma2_check > 2.0 * sigma2_vmc || sigma2_check < 0.5 * sigma2_vmc)
      throw std::runtime_error(
          "Safeguard failure: checkConfigurations variance out of [0.5, 2.0] * reference! Please report this bug.\n");
    app_log() << "  Execution time = " << std::setprecision(4) << t2.elapsed() << std::endl;
  }
  app_log() << "  </log>" << std::endl;
  app_log() << "</opt>" << std::endl;
  app_log() << "<opt stage=\"main\" walkers=\"" << optTarget->getNumSamples() << "\">" << std::endl;
  app_log() << "  <log>" << std::endl;
  t1.restart();
}

#ifdef HAVE_LMY_ENGINE
void QMCLinearOptimize::engine_start(cqmc::engine::LMYEngine<ValueType>* EngineObj,
                                     DescentEngine& descentEngineObj,
                                     std::string MinMethod)
{
  app_log() << "entering engine_start function" << std::endl;

  {
    //generate samples
    ScopedTimer local(generate_samples_timer_);
    generateSamples();
    //store active number of walkers
    NumOfVMCWalkers = W.getActiveWalkers();
  }

  app_log() << "<opt stage=\"setup\">" << std::endl;
  app_log() << "  <log>" << std::endl;

  // reset the root name
  optTarget->setRootName(RootName);
  optTarget->setWaveFunctionNode(wfNode);
  app_log() << "     Reading configurations from h5FileRoot " << h5FileRoot << std::endl;
  {
    // get configuration from the previous run
    ScopedTimer local(initialize_timer_);
    Timer t2;
    optTarget->getConfigurations(h5FileRoot);
    optTarget->setRng(vmcEngine->getRngRefs());
    optTarget->engine_checkConfigurations(EngineObj, descentEngineObj,
                                          MinMethod); // computes derivative ratios and pass into engine
    app_log() << "  Execution time = " << std::setprecision(4) << t2.elapsed() << std::endl;
  }
  app_log() << "  </log>" << std::endl;
  app_log() << "</opt>" << std::endl;
  app_log() << "<opt stage=\"main\" walkers=\"" << optTarget->getNumSamples() << "\">" << std::endl;
  app_log() << "  <log>" << std::endl;
  t1.restart();
}
#endif

void QMCLinearOptimize::finish()
{
  MyCounter++;
  app_log() << "  Execution time = " << std::setprecision(4) << t1.elapsed() << std::endl;
  app_log() << "  </log>" << std::endl;

  if (optTarget->reportH5)
    optTarget->reportParametersH5();
  optTarget->reportParameters();


  int nw_removed = W.getActiveWalkers() - NumOfVMCWalkers;
  app_log() << "   Restore the number of walkers to " << NumOfVMCWalkers << ", removing " << nw_removed << " walkers."
            << std::endl;
  if (nw_removed > 0)
    W.destroyWalkers(nw_removed);
  else
    W.createWalkers(-nw_removed);
  app_log() << "</opt>" << std::endl;
  app_log() << "</optimization-report>" << std::endl;
}

void QMCLinearOptimize::generateSamples()
{
  app_log() << "<optimization-report>" << std::endl;
  vmcEngine->qmc_driver_mode.set(QMC_WARMUP, 1);
  //  vmcEngine->run();
  //  vmcEngine->setValue("blocks",nBlocks);
  //  app_log() << "  Execution time = " << std::setprecision(4) << t1.elapsed() << std::endl;
  //  app_log() << "</vmc>" << std::endl;
  //}
  //     if (W.getActiveWalkers()>NumOfVMCWalkers)
  //     {
  //         W.destroyWalkers(W.getActiveWalkers()-NumOfVMCWalkers);
  //         app_log() << "  QMCLinearOptimize::generateSamples removed walkers." << std::endl;
  //         app_log() << "  Number of Walkers per node " << W.getActiveWalkers() << std::endl;
  //     }
  vmcEngine->qmc_driver_mode.set(QMC_OPTIMIZE, 1);
  vmcEngine->qmc_driver_mode.set(QMC_WARMUP, 0);
  //vmcEngine->setValue("recordWalkers",1);//set record
  vmcEngine->setValue("current", 0); //reset CurrentStep
  app_log() << "<vmc stage=\"main\" blocks=\"" << nBlocks << "\">" << std::endl;
  t1.restart();
  //     W.reset();
  branchEngine->flush(0);
  branchEngine->reset();
  vmcEngine->run();
  app_log() << "  Execution time = " << std::setprecision(4) << t1.elapsed() << std::endl;
  app_log() << "</vmc>" << std::endl;
  h5FileRoot = RootName;
}

QMCLinearOptimize::RealType QMCLinearOptimize::getLowestEigenvector(Matrix<RealType>& A,
                                                                    Matrix<RealType>& B,
                                                                    std::vector<RealType>& ev)
{
  int Nl(ev.size());
  //Tested the single eigenvalue speed and It was no faster.
  //segfault issues with single eigenvalue problem for some machines
  //  bool singleEV(false);
  //  if (singleEV)
  //  {
  /*
    Matrix<double> TAU(Nl,Nl);
    int INFO;
    int LWORK(-1);
    std::vector<RealType> WORK(1);
    //optimal work size
    dgeqrf( &Nl, &Nl, B.data(), &Nl, TAU.data(), &WORK[0], &LWORK, &INFO);
    LWORK=int(WORK[0]);
    WORK.resize(LWORK);
    //QR factorization of S, or H2 matrix. to be applied to H before solve.
    dgeqrf( &Nl, &Nl, B.data(), &Nl, TAU.data(), &WORK[0], &LWORK, &INFO);
    char SIDE('L');
    char TRANS('T');
    LWORK=-1;
    //optimal work size
    dormqr(&SIDE, &TRANS, &Nl, &Nl, &Nl, B.data(), &Nl, TAU.data(), A.data(), &Nl, &WORK[0], &LWORK, &INFO);
    LWORK=int(WORK[0]);
    WORK.resize(LWORK);
    //Apply Q^T to H
    dormqr(&SIDE, &TRANS, &Nl, &Nl, &Nl, B.data(), &Nl, TAU.data(), A.data(), &Nl, &WORK[0], &LWORK, &INFO);
    //now we have a pair (A,B)=(Q^T*H,Q^T*S) where B is upper triangular and A is general matrix.
    //reduce the matrix pair to generalized upper Hesenberg form
    char COMPQ('N'), COMPZ('I');
    int ILO(1);
    int LDQ(Nl);
    Matrix<double> Z(Nl,Nl), Q(Nl,LDQ); //starts as unit matrix
    for (int zi=0; zi<Nl; zi++)
      Z(zi,zi)=1;
    dgghrd(&COMPQ, &COMPZ, &Nl, &ILO, &Nl, A.data(), &Nl, B.data(), &Nl, Q.data(), &LDQ, Z.data(), &Nl, &INFO);
    //Take the pair and reduce to shur form and get eigenvalues
    std::vector<RealType> alphar(Nl),alphai(Nl),beta(Nl);
    char JOB('S');
    COMPQ='N';
    COMPZ='V';
    LWORK=-1;
    //get optimal work size
    dhgeqz(&JOB, &COMPQ, &COMPZ, &Nl, &ILO, &Nl, A.data(), &Nl, B.data(), &Nl, &alphar[0], &alphai[0], &beta[0], Q.data(), &LDQ, Z.data(), &Nl, &WORK[0], &LWORK, &INFO);
    LWORK=int(WORK[0]);
    WORK.resize(LWORK);
    dhgeqz(&JOB, &COMPQ, &COMPZ, &Nl, &ILO, &Nl, A.data(), &Nl, B.data(), &Nl, &alphar[0], &alphai[0], &beta[0], Q.data(), &LDQ, Z.data(), &Nl, &WORK[0], &LWORK, &INFO);
    //find the best eigenvalue
    std::vector<std::pair<RealType,int> > mappedEigenvalues(Nl);
    for (int i=0; i<Nl; i++)
    {
      RealType evi(alphar[i]/beta[i]);
      if (std::abs(evi)<1e10)
      {
        mappedEigenvalues[i].first=evi;
        mappedEigenvalues[i].second=i;
      }
      else
      {
        mappedEigenvalues[i].first=1e100;
        mappedEigenvalues[i].second=i;
      }
    }
    std::sort(mappedEigenvalues.begin(),mappedEigenvalues.end());
    int BestEV(mappedEigenvalues[0].second);
//                   now we rearrange the  the matrices
    if (BestEV!=0)
    {
      bool WANTQ(false);
      bool WANTZ(true);
      int ILST(1);
      int IFST(BestEV+1);
      LWORK=-1;
      dtgexc(&WANTQ, &WANTZ, &Nl, A.data(), &Nl, B.data(), &Nl, Q.data(), &Nl, Z.data(), &Nl, &IFST, &ILST, &WORK[0], &LWORK, &INFO);
      LWORK=int(WORK[0]);
      WORK.resize(LWORK);
      dtgexc(&WANTQ, &WANTZ, &Nl, A.data(), &Nl, B.data(), &Nl, Q.data(), &Nl, Z.data(), &Nl, &IFST, &ILST, &WORK[0], &LWORK, &INFO);
    }
    //now we compute the eigenvector
    SIDE='R';
    char HOWMNY('S');
    int M(0);
    Matrix<double> Z_I(Nl,Nl);
    bool SELECT[Nl];
    for (int zi=0; zi<Nl; zi++)
      SELECT[zi]=false;
    SELECT[0]=true;
    WORK.resize(6*Nl);
    dtgevc(&SIDE, &HOWMNY, &SELECT[0], &Nl, A.data(), &Nl, B.data(), &Nl, Q.data(), &LDQ, Z_I.data(), &Nl, &Nl, &M, &WORK[0], &INFO);
    std::vector<RealType> evec(Nl,0);
    for (int i=0; i<Nl; i++)
      for (int j=0; j<Nl; j++)
        evec[i] += Z(j,i)*Z_I(0,j);
    for (int i=0; i<Nl; i++)
      ev[i] = evec[i]/evec[0];
//     for (int i=0; i<Nl; i++) app_log()<<ev[i]<<" ";
//     app_log()<< std::endl;
    return mappedEigenvalues[0].first;
    */
  // a fake return to reduce warning.
  //    return RealType(0.0);
  //  }
  //  else
  //  {
  // OLD ROUTINE. CALCULATES ALL EIGENVECTORS
  //   Getting the optimal worksize
  char jl('N');
  char jr('V');
  std::vector<RealType> alphar(Nl), alphai(Nl), beta(Nl);
  Matrix<RealType> eigenT(Nl, Nl);
  int info;
  int lwork(-1);
  std::vector<RealType> work(1);
  RealType tt(0);
  int t(1);
  LAPACK::ggev(&jl, &jr, &Nl, A.data(), &Nl, B.data(), &Nl, &alphar[0], &alphai[0], &beta[0], &tt, &t, eigenT.data(),
               &Nl, &work[0], &lwork, &info);
  lwork = int(work[0]);
  work.resize(lwork);
  //~ //Get an estimate of E_lin
  //~ Matrix<RealType> H_tmp(HamT);
  //~ Matrix<RealType> S_tmp(ST);
  //~ dggev(&jl, &jr, &Nl, H_tmp.data(), &Nl, S_tmp.data(), &Nl, &alphar[0], &alphai[0], &beta[0],&tt,&t, eigenT.data(), &Nl, &work[0], &lwork, &info);
  //~ RealType E_lin(alphar[0]/beta[0]);
  //~ int e_min_indx(0);
  //~ for (int i=1; i<Nl; i++)
  //~ if (E_lin>(alphar[i]/beta[i]))
  //~ {
  //~ E_lin=alphar[i]/beta[i];
  //~ e_min_indx=i;
  //~ }
  LAPACK::ggev(&jl, &jr, &Nl, A.data(), &Nl, B.data(), &Nl, &alphar[0], &alphai[0], &beta[0], &tt, &t, eigenT.data(),
               &Nl, &work[0], &lwork, &info);
  if (info != 0)
  {
    APP_ABORT("Invalid Matrix Diagonalization Function!");
  }
  std::vector<std::pair<RealType, int>> mappedEigenvalues(Nl);
  for (int i = 0; i < Nl; i++)
  {
    RealType evi(alphar[i] / beta[i]);
    if (std::abs(evi) < 1e10)
    {
      mappedEigenvalues[i].first  = evi;
      mappedEigenvalues[i].second = i;
    }
    else
    {
      mappedEigenvalues[i].first  = std::numeric_limits<RealType>::max();
      mappedEigenvalues[i].second = i;
    }
  }
  std::sort(mappedEigenvalues.begin(), mappedEigenvalues.end());
  for (int i = 0; i < Nl; i++)
    ev[i] = eigenT(mappedEigenvalues[0].second, i) / eigenT(mappedEigenvalues[0].second, 0);
  return mappedEigenvalues[0].first;
  //  }
}


QMCLinearOptimize::RealType QMCLinearOptimize::getLowestEigenvector(Matrix<RealType>& A, std::vector<RealType>& ev)
{
  int Nl(ev.size());
  //Tested the single eigenvalue speed and It was no faster.
  //segfault issues with single eigenvalue problem for some machines
  //  bool singleEV(false);
  //     if (singleEV)
  //     {
  //         Matrix<double> TAU(Nl,Nl);
  //         int INFO;
  //         int LWORK(-1);
  //         std::vector<RealType> WORK(1);
  //         //optimal work size
  //         dgeqrf( &Nl, &Nl, B.data(), &Nl, TAU.data(), &WORK[0], &LWORK, &INFO);
  //         LWORK=int(WORK[0]);
  //         WORK.resize(LWORK);
  //         //QR factorization of S, or H2 matrix. to be applied to H before solve.
  //         dgeqrf( &Nl, &Nl, B.data(), &Nl, TAU.data(), &WORK[0], &LWORK, &INFO);
  //
  //         char SIDE('L');
  //         char TRANS('T');
  //         LWORK=-1;
  //         //optimal work size
  //         dormqr(&SIDE, &TRANS, &Nl, &Nl, &Nl, B.data(), &Nl, TAU.data(), A.data(), &Nl, &WORK[0], &LWORK, &INFO);
  //         LWORK=int(WORK[0]);
  //         WORK.resize(LWORK);
  //         //Apply Q^T to H
  //         dormqr(&SIDE, &TRANS, &Nl, &Nl, &Nl, B.data(), &Nl, TAU.data(), A.data(), &Nl, &WORK[0], &LWORK, &INFO);
  //
  //         //now we have a pair (A,B)=(Q^T*H,Q^T*S) where B is upper triangular and A is general matrix.
  //         //reduce the matrix pair to generalized upper Hesenberg form
  //         char COMPQ('N'), COMPZ('I');
  //         int ILO(1);
  //         int LDQ(Nl);
  //         Matrix<double> Z(Nl,Nl), Q(Nl,LDQ); //starts as unit matrix
  //         for (int zi=0; zi<Nl; zi++) Z(zi,zi)=1;
  //         dgghrd(&COMPQ, &COMPZ, &Nl, &ILO, &Nl, A.data(), &Nl, B.data(), &Nl, Q.data(), &LDQ, Z.data(), &Nl, &INFO);
  //
  //         //Take the pair and reduce to shur form and get eigenvalues
  //         std::vector<RealType> alphar(Nl),alphai(Nl),beta(Nl);
  //         char JOB('S');
  //         COMPQ='N';
  //         COMPZ='V';
  //         LWORK=-1;
  //         //get optimal work size
  //         dhgeqz(&JOB, &COMPQ, &COMPZ, &Nl, &ILO, &Nl, A.data(), &Nl, B.data(), &Nl, &alphar[0], &alphai[0], &beta[0], Q.data(), &LDQ, Z.data(), &Nl, &WORK[0], &LWORK, &INFO);
  //         LWORK=int(WORK[0]);
  //         WORK.resize(LWORK);
  //         dhgeqz(&JOB, &COMPQ, &COMPZ, &Nl, &ILO, &Nl, A.data(), &Nl, B.data(), &Nl, &alphar[0], &alphai[0], &beta[0], Q.data(), &LDQ, Z.data(), &Nl, &WORK[0], &LWORK, &INFO);
  //         //find the best eigenvalue
  //         std::vector<std::pair<RealType,int> > mappedEigenvalues(Nl);
  //         for (int i=0; i<Nl; i++)
  //         {
  //             RealType evi(alphar[i]/beta[i]);
  //             if (std::abs(evi)<1e10)
  //             {
  //                 mappedEigenvalues[i].first=evi;
  //                 mappedEigenvalues[i].second=i;
  //             }
  //             else
  //             {
  //                 mappedEigenvalues[i].first=1e100;
  //                 mappedEigenvalues[i].second=i;
  //             }
  //         }
  //         std::sort(mappedEigenvalues.begin(),mappedEigenvalues.end());
  //         int BestEV(mappedEigenvalues[0].second);
  //
  // //                   now we rearrange the  the matrices
  //         if (BestEV!=0)
  //         {
  //             bool WANTQ(false);
  //             bool WANTZ(true);
  //             int ILST(1);
  //             int IFST(BestEV+1);
  //             LWORK=-1;
  //
  //             dtgexc(&WANTQ, &WANTZ, &Nl, A.data(), &Nl, B.data(), &Nl, Q.data(), &Nl, Z.data(), &Nl, &IFST, &ILST, &WORK[0], &LWORK, &INFO);
  //             LWORK=int(WORK[0]);
  //             WORK.resize(LWORK);
  //             dtgexc(&WANTQ, &WANTZ, &Nl, A.data(), &Nl, B.data(), &Nl, Q.data(), &Nl, Z.data(), &Nl, &IFST, &ILST, &WORK[0], &LWORK, &INFO);
  //         }
  //         //now we compute the eigenvector
  //         SIDE='R';
  //         char HOWMNY('S');
  //         int M(0);
  //         Matrix<double> Z_I(Nl,Nl);
  //         bool SELECT[Nl];
  //         for (int zi=0; zi<Nl; zi++) SELECT[zi]=false;
  //         SELECT[0]=true;
  //
  //         WORK.resize(6*Nl);
  //         dtgevc(&SIDE, &HOWMNY, &SELECT[0], &Nl, A.data(), &Nl, B.data(), &Nl, Q.data(), &LDQ, Z_I.data(), &Nl, &Nl, &M, &WORK[0], &INFO);
  //
  //         std::vector<RealType> evec(Nl,0);
  //         for (int i=0; i<Nl; i++) for (int j=0; j<Nl; j++) evec[i] += Z(j,i)*Z_I(0,j);
  //         for (int i=0; i<Nl; i++) ev[i] = evec[i]/evec[0];
  // //     for (int i=0; i<Nl; i++) app_log()<<ev[i]<<" ";
  // //     app_log()<< std::endl;
  //         return mappedEigenvalues[0].first;
  //     }
  //     else
  //     {
  // // OLD ROUTINE. CALCULATES ALL EIGENVECTORS
  // //   Getting the optimal worksize
  RealType zerozero = A(0, 0);
  char jl('N');
  char jr('V');
  std::vector<RealType> alphar(Nl), alphai(Nl), beta(Nl);
  Matrix<RealType> eigenT(Nl, Nl);
  Matrix<RealType> eigenD(Nl, Nl);
  int info;
  int lwork(-1);
  std::vector<RealType> work(1);
  LAPACK::geev(&jl, &jr, &Nl, A.data(), &Nl, &alphar[0], &alphai[0], eigenD.data(), &Nl, eigenT.data(), &Nl, &work[0],
               &lwork, &info);
  lwork = int(work[0]);
  work.resize(lwork);
  //~ //Get an estimate of E_lin
  //~ Matrix<RealType> H_tmp(HamT);
  //~ Matrix<RealType> S_tmp(ST);
  //~ dggev(&jl, &jr, &Nl, H_tmp.data(), &Nl, S_tmp.data(), &Nl, &alphar[0], &alphai[0], &beta[0],&tt,&t, eigenT.data(), &Nl, &work[0], &lwork, &info);
  //~ RealType E_lin(alphar[0]/beta[0]);
  //~ int e_min_indx(0);
  //~ for (int i=1; i<Nl; i++)
  //~ if (E_lin>(alphar[i]/beta[i]))
  //~ {
  //~ E_lin=alphar[i]/beta[i];
  //~ e_min_indx=i;
  //~ }
  LAPACK::geev(&jl, &jr, &Nl, A.data(), &Nl, &alphar[0], &alphai[0], eigenD.data(), &Nl, eigenT.data(), &Nl, &work[0],
               &lwork, &info);
  if (info != 0)
  {
    APP_ABORT("Invalid Matrix Diagonalization Function!");
  }
  std::vector<std::pair<RealType, int>> mappedEigenvalues(Nl);
  for (int i = 0; i < Nl; i++)
  {
    RealType evi(alphar[i]);
    if ((evi < zerozero) && (evi > (zerozero - 1e2)))
    {
      mappedEigenvalues[i].first  = (evi - zerozero + 2.0) * (evi - zerozero + 2.0);
      mappedEigenvalues[i].second = i;
    }
    else
    {
      mappedEigenvalues[i].first  = std::numeric_limits<RealType>::max();
      mappedEigenvalues[i].second = i;
    }
  }
  std::sort(mappedEigenvalues.begin(), mappedEigenvalues.end());
  //         for (int i=0; i<4; i++) app_log()<<i<<": "<<alphar[mappedEigenvalues[i].second]<< std::endl;
  for (int i = 0; i < Nl; i++)
    ev[i] = eigenT(mappedEigenvalues[0].second, i) / eigenT(mappedEigenvalues[0].second, 0);
  return alphar[mappedEigenvalues[0].second];
  //     }
}
bool QMCLinearOptimize::nonLinearRescale(std::vector<RealType>& dP, Matrix<RealType>& S)
{
  RealType rescale = getNonLinearRescale(dP, S);
  for (int i = 1; i < dP.size(); i++)
    dP[i] *= rescale;
  return true;
}


void QMCLinearOptimize::getNonLinearRange(int& first, int& last)
{
  std::vector<int> types;
  optTarget->getParameterTypes(types);
  first = 0;
  last  = types.size();
  //assume all non-linear coeffs are together.
  if (types[0] == optimize::LINEAR_P)
  {
    int i(0);
    while (i < types.size())
    {
      if (types[i] == optimize::LINEAR_P)
        first = i;
      i++;
    }
    first++;
  }
  else
  {
    int i(types.size() - 1);
    while (i >= 0)
    {
      if (types[i] == optimize::LINEAR_P)
        last = i;
      i--;
    }
  }
  //     returns the number of non-linear parameters.
  //    app_log()<<"line params: "<<first<<" "<<last<< std::endl;
}

QMCLinearOptimize::RealType QMCLinearOptimize::getNonLinearRescale(std::vector<RealType>& dP, Matrix<RealType>& S)
{
  int first(0), last(0);
  getNonLinearRange(first, last);
  if (first == last)
    return 1.0;
  RealType rescale(1.0);
  RealType xi(0.5);
  RealType D(0.0);
  for (int i = first; i < last; i++)
    for (int j = first; j < last; j++)
      D += S(i + 1, j + 1) * dP[i + 1] * dP[j + 1];
  rescale = (1 - xi) * D / ((1 - xi) + xi * std::sqrt(1 + D));
  rescale = 1.0 / (1.0 - rescale);
  //     app_log()<<"rescale: "<<rescale<< std::endl;
  return rescale;
}

void QMCLinearOptimize::orthoScale(std::vector<RealType>& dP, Matrix<RealType>& S)
{
  //     int first(0),last(0);
  //     getNonLinearRange(first,last);
  //     if (first==last) return;
  int x(dP.size());
  Matrix<RealType> T(S);
  std::vector<RealType> nP(dP);
  Matrix<RealType> lS(x, x);
  for (int i = 0; i < x; i++)
    for (int j = 0; j < x; j++)
      lS(i, j) = S(i + 1, j + 1);
  RealType Det = invert_matrix(lS, true);
  for (int i = 0; i < x; i++)
  {
    dP[i] = 0;
    for (int j = 0; j < x; j++)
    {
      dP[i] += nP[j] * lS(i, j);
    }
  }
  RealType rs = getNonLinearRescale(dP, T);
  for (int i = 0; i < x; i++)
    dP[i] *= rs;
  for (int i = 0; i < dP.size(); i++)
    app_log() << dP[i] << " ";
  app_log() << std::endl;
  //     RealType D(0.0);
  //     for (int i=first; i<last; i++) D += 2.0*S(i+1,0)*dP[i];
  //     for (int i=first; i<last; i++) for (int j=first; j<last; j++) D += S(i+1,j+1)*dP[i]*dP[j];
  //     app_log()<<D<< std::endl;
  //
  //
  //     RealType rescale = 0.5*D/(0.5 + 0.5*std::sqrt(1.0+D));
  //     rescale = 1.0/(1.0-rescale);
  //     app_log()<<rescale<< std::endl;
  //     for (int i=0; i<dP.size(); i++) dP[i] *= rescale;
}

/** Parses the xml input file for parameter definitions for the wavefunction optimization.
* @param q current xmlNode
* @return true if successful
*/
bool QMCLinearOptimize::put(xmlNodePtr q)
{
  std::string useGPU("no");
  std::string vmcMove("pbyp");
  OhmmsAttributeSet oAttrib;
  oAttrib.add(useGPU, "gpu");
  oAttrib.add(vmcMove, "move");
  oAttrib.put(q);
  optNode        = q;
  xmlNodePtr cur = optNode->children;
  int pid        = OHMMS::Controller->rank();
  while (cur != NULL)
  {
    std::string cname((const char*)(cur->name));
    if (cname == "mcwalkerset")
    {
      mcwalkerNodePtr.push_back(cur);
    }
    cur = cur->next;
  }
  //no walkers exist, add 10
  if (W.getActiveWalkers() == 0)
    addWalkers(omp_get_max_threads());
  NumOfVMCWalkers = W.getActiveWalkers();
  bool success    = true;
  //allways reset optTarget
#if defined(QMC_CUDA)
  if (useGPU == "yes")
    optTarget = std::make_unique<QMCCostFunctionCUDA>(W, Psi, H, myComm);
  else
#endif
    optTarget = std::make_unique<QMCCostFunction>(W, Psi, H, myComm);
  optTarget->setStream(&app_log());
  success = optTarget->put(q);

  //create VMC engine
  if (vmcEngine == 0)
  {
#if defined(QMC_CUDA)
    if (useGPU == "yes")
      vmcEngine = std::make_unique<VMCcuda>(W, Psi, H, myComm, false);
    else
#endif
      vmcEngine = std::make_unique<VMC>(W, Psi, H, myComm, false);
    vmcEngine->setUpdateMode(vmcMove[0] == 'p');
  }

  vmcEngine->setStatus(RootName, h5FileRoot, AppendRun);
  vmcEngine->process(optNode);
  return success;
}

bool QMCLinearOptimize::fitMappedStabilizers(std::vector<std::pair<RealType, RealType>>& mappedStabilizers,
                                             RealType& XS,
                                             RealType& val,
                                             RealType tooBig)
{
  int nms(0);
  for (int i = 0; i < mappedStabilizers.size(); i++)
    if (mappedStabilizers[i].second == mappedStabilizers[i].second)
      nms++;
  bool SuccessfulFit(false);
  if (nms >= 5)
  {
    //Quartic fit the stabilizers we have tried and try to choose the best we can
    std::vector<RealType> Y(nms), Coefs(5);
    Matrix<RealType> X(nms, 5);
    for (int i = 0; i < nms; i++)
      if (mappedStabilizers[i].second == mappedStabilizers[i].second)
      {
        X(i, 0) = 1.0;
        X(i, 1) = mappedStabilizers[i].first;
        X(i, 2) = std::pow(mappedStabilizers[i].first, 2);
        X(i, 3) = std::pow(mappedStabilizers[i].first, 3);
        X(i, 4) = std::pow(mappedStabilizers[i].first, 4);
        Y[i]    = mappedStabilizers[i].second;
      }
    LinearFit(Y, X, Coefs);
    RealType Xmin = QuarticMinimum(Coefs);
    val           = 0;
    for (int i = 0; i < 5; i++)
      val += std::pow(Xmin, i) * Coefs[i];
    app_log() << "quartic Fit min: " << Xmin << " val: " << val << std::endl;
    ;
    //         for (int i=0; i<5; i++) app_log()<<Coefs[i]<<" ";
    //         app_log()<< std::endl;
    SuccessfulFit = true;
    for (int i = 0; i < nms; i++)
      if (mappedStabilizers[i].second == mappedStabilizers[i].second)
        if (val > mappedStabilizers[i].second)
          SuccessfulFit = false;
    if (Xmin > tooBig)
      SuccessfulFit = false;
    if (SuccessfulFit)
      XS = Xmin;
  }
  else if (nms >= 3)
  {
    //Quadratic fit the stabilizers we have tried and try to choose the best we can
    std::sort(mappedStabilizers.begin(), mappedStabilizers.end());
    std::vector<RealType> Y(nms), Coefs(3);
    Matrix<RealType> X(nms, 3);
    for (int i = 0; i < nms; i++)
      if (mappedStabilizers[i].second == mappedStabilizers[i].second)
      {
        X(i, 0) = 1.0;
        X(i, 1) = mappedStabilizers[i].first;
        X(i, 2) = std::pow(mappedStabilizers[i].first, 2);
        Y[i]    = mappedStabilizers[i].second;
      }
    LinearFit(Y, X, Coefs);
    //extremum really.
    RealType Xmin = -0.5 * Coefs[1] / Coefs[2];
    val           = 0;
    for (int i = 0; i < 3; i++)
      val += std::pow(Xmin, i) * Coefs[i];
    app_log() << "quadratic Fit min: " << Xmin << " val: " << val << std::endl;
    //         for (int i=0; i<3; i++) app_log()<<Coefs[i]<<" ";
    //         app_log()<< std::endl;
    SuccessfulFit = true;
    if (Xmin > tooBig)
      SuccessfulFit = false;
    for (int i = 0; i < nms; i++)
      if (mappedStabilizers[i].second == mappedStabilizers[i].second)
        if (val > mappedStabilizers[i].second)
          SuccessfulFit = false;
    if (SuccessfulFit)
      XS = Xmin;
  }
  return SuccessfulFit;
}
} // namespace qmcplusplus
