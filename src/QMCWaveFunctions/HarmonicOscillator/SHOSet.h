//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jaron T. Krogel, krogeljt@ornl.gov, Oak Ridge National Laboratory
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jaron T. Krogel, krogeljt@ornl.gov, Oak Ridge National Laboratory
//////////////////////////////////////////////////////////////////////////////////////


#ifndef QMCPLUSPLUS_SHOSET_H
#define QMCPLUSPLUS_SHOSET_H

#include "QMCWaveFunctions/SPOSet.h"
#include "QMCWaveFunctions/SPOInfo.h"


namespace qmcplusplus
{
struct SHOState : public SPOInfo
{
  TinyVector<int, DIM> quantum_number;

  SHOState()
  {
    quantum_number = -1;
    energy         = 0.0;
  }

  ~SHOState() override {}

  inline void set(TinyVector<int, DIM> qn, RealType e)
  {
    quantum_number = qn;
    energy         = e;
  }

  inline void sho_report(const std::string& pad = "") const
  {
    app_log() << pad << "qn=" << quantum_number << "  e=" << energy << std::endl;
  }
};


struct SHOSet : public SPOSet
{
  using value_type = ValueMatrix::value_type;
  using grad_type  = GradMatrix::value_type;

  RealType length;
  PosType center;

  int nmax;
  TinyVector<int, DIM> qn_max;
  std::vector<SHOState> state_info;
  std::vector<RealType> prefactors;
  Array<RealType, 2> hermite;
  Array<RealType, 2> bvalues;
  Array<RealType, 2> d0_values;
  Array<RealType, 2> d1_values;
  Array<RealType, 2> d2_values;

  //construction/destruction
  SHOSet(const std::string& my_name, RealType l, PosType c, const std::vector<SHOState*>& sho_states);

  ~SHOSet() override;

  std::string getClassName() const override { return "SHOSet"; }

  void initialize();

  //SPOSet interface methods
  std::unique_ptr<SPOSet> makeClone() const override;

  void evaluateValue(const ParticleSet& P, int iat, ValueVector& psi) override;

  void evaluateVGL(const ParticleSet& P, int iat, ValueVector& psi, GradVector& dpsi, ValueVector& d2psi) override;

  void evaluate_notranspose(const ParticleSet& P,
                            int first,
                            int last,
                            ValueMatrix& logdet,
                            GradMatrix& dlogdet,
                            ValueMatrix& d2logdet) override;


  //local functions
  void evaluate_v(PosType r, ValueVector& psi);
  void evaluate_vgl(PosType r, ValueVector& psi, GradVector& dpsi, ValueVector& d2psi);
  void evaluate_hermite(const PosType& xpos);
  void evaluate_d0(const PosType& xpos, ValueVector& psi);
  void evaluate_d1(const PosType& xpos, ValueVector& psi, GradVector& dpsi);
  void evaluate_d2(const PosType& xpos, ValueVector& psi, ValueVector& d2psi);
  void report(const std::string& pad = "") const override;
  void test_derivatives();
  void test_overlap();
  void evaluate_check(PosType r, ValueVector& psi, GradVector& dpsi, ValueVector& d2psi);

  //empty methods
  /// number of orbitals is determined only by initial request
  inline void setOrbitalSetSize(int norbs) override {}

  ///unimplemented functions call this to abort
  inline void not_implemented(const std::string& method)
  {
    APP_ABORT("SHOSet::" + method + " has not been implemented.");
  }


  //methods to be implemented in the future (possibly)
  void evaluateThirdDeriv(const ParticleSet& P, int first, int last, GGGMatrix& dddlogdet) override;
  void evaluate_notranspose(const ParticleSet& P,
                            int first,
                            int last,
                            ValueMatrix& logdet,
                            GradMatrix& dlogdet,
                            HessMatrix& ddlogdet) override;
  void evaluate_notranspose(const ParticleSet& P,
                            int first,
                            int last,
                            ValueMatrix& logdet,
                            GradMatrix& dlogdet,
                            HessMatrix& ddlogdet,
                            GGGMatrix& dddlogdet) override;
  void evaluateGradSource(const ParticleSet& P,
                          int first,
                          int last,
                          const ParticleSet& source,
                          int iat_src,
                          GradMatrix& gradphi) override;
  void evaluateGradSource(const ParticleSet& P,
                          int first,
                          int last,
                          const ParticleSet& source,
                          int iat_src,
                          GradMatrix& dphi,
                          HessMatrix& ddphi,
                          GradMatrix& dlapl_phi) override;
};

} // namespace qmcplusplus


#endif
