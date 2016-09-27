//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Miguel Morales, moralessilva2@llnl.gov, Lawrence Livermore National Laboratory
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////
    
    

#ifndef QMCPLUSPLUS_EINSPLINE_SET_LOCAL_H
#define QMCPLUSPLUS_EINSPLINE_SET_LOCAL_H

#include "QMCWaveFunctions/EinsplineOrb.h"
#include "QMCWaveFunctions/EinsplineSet.h"

namespace qmcplusplus
{

class EinsplineSetLocal : public EinsplineSet
{
  friend class EinsplineSetBuilder;
protected:
  /////////////////////
  // Orbital storage //
  /////////////////////
  /// Store the orbital objects.  Using template class allows us to
  /// avoid making separate real and complex versions of this class.
  //std::vector<EinsplineOrb<ValueType,OHMMS_DIM>*> Orbitals;
  std::vector<EinsplineOrb<std::complex<double>,OHMMS_DIM>*> Orbitals;

public:
  SPOSetBase* makeClone() const;
  void evaluate(const ParticleSet& P, int iat, ValueVector_t& psi);
  void evaluate(const ParticleSet& P, int iat,
                ValueVector_t& psi, GradVector_t& dpsi, ValueVector_t& d2psi);
  void evaluate(const ParticleSet& P, int iat,
                ValueVector_t& psi, GradVector_t& dpsi, HessVector_t& gg_psi)
  {
    APP_ABORT("Need specialization of evaluate(iat) for HessVector. \n");
  }
  void evaluate_notranspose(const ParticleSet& P, int first, int last,
                            ValueMatrix_t& psi, GradMatrix_t& dpsi,
                            ValueMatrix_t& d2psi);
  void evaluate_notranspose(const ParticleSet& P, int first, int last,
                            ValueMatrix_t& logdet, GradMatrix_t& dlogdet,
                            HessMatrix_t& grad_grad_logdet);
  void evaluate_notranspose(const ParticleSet& P, int first, int last,
                            ValueMatrix_t& logdet, GradMatrix_t& dlogdet,
                            HessMatrix_t& grad_grad_logdet,
                            GGGMatrix_t& grad_grad_grad_logdet);

  void resetParameters(const opt_variables_type& active);

  EinsplineSetLocal()
  {
    className = "EinsplineSetLocal";
  }
};
}
#endif
