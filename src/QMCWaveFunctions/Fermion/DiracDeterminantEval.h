//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2018 QMCPACK developers.
//
// File developed by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Laboratory
//   partially refactored from DiracDeterminantBase.h
//
// File created by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Laboratory
//////////////////////////////////////////////////////////////////////////////////////

#ifndef QMCPLUSPLUS_DIRAC_DETERMINANT_EVAL_H
#define QMCPLUSPLUS_DIRAC_DETERMINANT_EVAL_H

#include "Configuration.h"
#include "QMCWaveFunctions/Batching.h"
#include "QMCWaveFunctions/SPOSetSingle.h"
#include "QMCWaveFunctions/SPOSetBatched.h"
#include "QMCWaveFunctions/Fermion/determinant_update.h"
#include "Particle/MCWalkerConfiguration.h"
#include "QMCWaveFunctions/SPOSetTypeAliases.h"
#include "QMCWaveFunctions/NLjob.h"
#include "QMCWaveFunctions/WaveFunctionComponentTypeAliases.h"
namespace qmcplusplus
{

class DiracDeterminantEvalDefault
{
public:
  using QMCT = QMCTraits;
  using SSTA = SPOSetTypeAliases;
  using WFCA = WaveFunctionComponentTypeAliases;
  using Walker_t =  WFCA::Walker_t;

  
  ///total number of particles
  int NP;
  ///number of single-particle orbitals which belong to this Dirac determinant
  int NumOrbitals;
  ///number of particles which belong to this Dirac determinant
  int NumPtcls;
  ///index of the first particle with respect to the particle set
  int FirstIndex;
  ///index of the last particle with respect to the particle set
  int LastIndex;
  ///index of the particle (or row)
  int WorkingIndex;
  ///a set of single-particle orbitals used to fill in the  values of the matrix
  //SPOSetPtr Phi;

  /////Current determinant value
  //ValueType CurrentDet;
  /// psiM(j,i) \f$= \psi_j({\bf r}_i)\f$
  ValueMatrix_t psiM, psiM_temp;

  /// memory pool for temporal data
  aligned_vector<ValueType> memoryPool;

  /// temporary container for testing
  ValueMatrix_t psiMinv;

  /// dpsiM(i,j) \f$= \nabla_i \psi_j({\bf r}_i)\f$
  GradMatrix_t  dpsiM, dpsiM_temp;

  /// d2psiM(i,j) \f$= \nabla_i^2 \psi_j({\bf r}_i)\f$
  ValueMatrix_t d2psiM, d2psiM_temp;

  /// Used for force computations
  GradMatrix_t grad_source_psiM, grad_lapl_source_psiM;
  HessMatrix_t  grad_grad_source_psiM;
  
  GradMatrix_t phi_alpha_Minv, grad_phi_Minv;
  ValueMatrix_t lapl_phi_Minv;
  HessMatrix_t grad_phi_alpha_Minv;

  /// value of single-particle orbital for particle-by-particle update
  ValueVector_t psiV;
  GradVector_t dpsiV;
  ValueVector_t d2psiV;
  ValueVector_t workV1, workV2;
  GradVector_t workG;

  
  void abortNoSpecialize()
  {
    APP_ABORT("DiracDeterminantEvalDefault methods should not be reached");
  }
  
  virtual void addLog (MCWalkerConfiguration &W, std::vector<QMCT::RealType> &logPsi)
  { abortNoSpecialize(); }
  virtual void addGradient(MCWalkerConfiguration &W, int iat,
			   std::vector<QMCT::GradType> &grad)
  { abortNoSpecialize(); }
  virtual void calcGradient(MCWalkerConfiguration &W, int iat,
			    std::vector<QMCT::GradType> &grad)
  { abortNoSpecialize(); }
  virtual void ratio (MCWalkerConfiguration &W, int iat,
		      std::vector<QMCT::ValueType> &psi_ratios)
  { abortNoSpecialize(); }

  virtual void ratio (MCWalkerConfiguration &W, int iat,
		      std::vector<QMCT::ValueType> &psi_ratios,
		      std::vector<QMCT::GradType>  &grad)
  { abortNoSpecialize(); }
  virtual void ratio (MCWalkerConfiguration &W, int iat,
		      std::vector<QMCT::ValueType> &psi_ratios,
		      std::vector<QMCT::GradType>  &grad,
		      std::vector<QMCT::ValueType> &lapl)
  { abortNoSpecialize(); }
  virtual void ratio (std::vector<Walker_t*> &walkers, std::vector<int> &iatList,
		      std::vector<QMCT::PosType> &rNew, std::vector<QMCT::ValueType> &psi_ratios,
		      std::vector<QMCT::GradType>  &grad, std::vector<QMCT::ValueType> &lapl)
  { abortNoSpecialize(); }
  virtual  void NLratios (MCWalkerConfiguration &W,  std::vector<NLjob> &jobList,
			  std::vector<QMCT::PosType> &quadPoints, std::vector<QMCT::ValueType> &psi_ratios)
  { abortNoSpecialize(); }
  virtual void calcRatio (MCWalkerConfiguration &W, int iat,
			  std::vector<QMCT::ValueType> &psi_ratios,	std::vector<QMCT::GradType>  &grad,
                  std::vector<QMCT::ValueType> &lapl)
  { abortNoSpecialize(); }
  virtual void addRatio (MCWalkerConfiguration &W, int iat,
                 std::vector<QMCT::ValueType> &psi_ratios,	std::vector<QMCT::GradType>  &grad,
                 std::vector<QMCT::ValueType> &lapl)
  { abortNoSpecialize(); }

  virtual void gradLapl (MCWalkerConfiguration &W, SSTA::GradMatrix_t &grads,
			 SSTA::ValueMatrix_t &lapl)
  { abortNoSpecialize(); }

  virtual QMCT::RealType updateBuffer(ParticleSet& P, WFCA::WFBufferType& buf, bool fromscratch=false)
  { abortNoSpecialize();
    return QMCT::RealType();
  }
  
  virtual void recompute(ParticleSet& P)
  {
    APP_ABORT("Need specialization of DiracDetermiantEval::recompute(ParticleSet& P).\n");
  }
  
  virtual void recompute(MCWalkerConfiguration &W, bool firstTime)
  {
    std::cerr << "Need specialization of DiracDetermiantEval::recompute.\n";
    abort();
  }

  virtual void update (std::vector<Walker_t*> &walkers, int iat)
  {
    APP_ABORT("Need specialization of DiracDetermiantEval::update.\n");
  }

  virtual void update (const std::vector<Walker_t*> &walkers, const std::vector<int> &iatList)
  {
    APP_ABORT("Need specialization of DiracDetermiantEval::update.\n");
  }


  
};

template<Batching batching>
class DiracDeterminantEval : public DiracDeterminantEvalDefault
{

};

template<>
class DiracDeterminantEval<Batching::SINGLE> : public DiracDeterminantEvalDefault
{
};

  
}

#endif
