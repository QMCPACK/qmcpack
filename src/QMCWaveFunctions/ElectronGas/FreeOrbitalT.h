//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2022 QMCPACK developers.
//
// File developed by: Miguel Morales, moralessilva2@llnl.gov, Lawrence Livermore National Laboratory
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Jaron T. Krogel, krogeljt@ornl.gov, Oak Ridge National Laboratory
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//                    Yubo "Paul" Yang, yubo.paul.yang@gmail.com, CCQ @ Flatiron
//                    William F Godoy, godoywf@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////

#ifndef QMCPLUSPLUS_FREE_ORBITALT_H
#define QMCPLUSPLUS_FREE_ORBITALT_H

#include "QMCWaveFunctions/SPOSetT.h"

namespace qmcplusplus
{
template<class T>
class FreeOrbitalT : public SPOSetT<T>
{
public:
  using ValueVector = typename SPOSetT<T>::ValueVector;
  using GradVector  = typename SPOSetT<T>::GradVector;
  using HessVector  = typename SPOSetT<T>::HessVector;
  using ValueMatrix = typename SPOSetT<T>::ValueMatrix;
  using GradMatrix  = typename SPOSetT<T>::GradMatrix;
  using HessMatrix  = typename SPOSetT<T>::HessMatrix;
  using GGGMatrix   = typename SPOSetT<T>::GGGMatrix;
  using RealType    = typename SPOSetT<T>::RealType;
  using PosType     = typename SPOSetT<T>::PosType;
  using ValueType   = typename SPOSetT<T>::ValueType;

  FreeOrbitalT(const std::string& my_name, const std::vector<PosType>& kpts_cart);
  ~FreeOrbitalT();

  inline std::string getClassName() const final { return "FreeOrbital"; }

  // phi[i][j] is phi_j(r_i), i.e. electron i in orbital j
  //  i \in [first, last)
  void evaluate_notranspose(const ParticleSetT<T>& P,
                            int first,
                            int last,
                            ValueMatrix& phi,
                            GradMatrix& dphi,
                            ValueMatrix& d2phi) final;

  // plug r_i into all orbitals
  void evaluateVGL(const ParticleSetT<T>& P, int i, ValueVector& pvec, GradVector& dpvec, ValueVector& d2pvec) final;
  void evaluateValue(const ParticleSetT<T>& P, int iat, ValueVector& pvec) final;

  // hessian matrix is needed by backflow
  void evaluate_notranspose(const ParticleSetT<T>& P,
                            int first,
                            int last,
                            ValueMatrix& phi,
                            GradMatrix& dphi,
                            HessMatrix& d2phi_mat) final;

  // derivative of hessian is needed to optimize backflow
  void evaluate_notranspose(const ParticleSetT<T>& P,
                            int first,
                            int last,
                            ValueMatrix& phi,
                            GradMatrix& dphi,
                            HessMatrix& d2phi_mat,
                            GGGMatrix& d3phi_mat) override;

  void report(const std::string& pad) const override;
  // ---- begin required overrides
  std::unique_ptr<SPOSetT<T>> makeClone() const final { return std::make_unique<FreeOrbitalT<T>>(*this); }
  void setOrbitalSetSize(int norbs) final { throw std::runtime_error("not implemented"); }
  // required overrides end ----
private:
  const std::vector<PosType> kvecs; // kvecs vectors
  const int mink;                   // minimum k index
  const int maxk;                   // maximum number of kvecs vectors
  std::vector<RealType> k2neg;      // minus kvecs^2
};

} // namespace qmcplusplus
#endif
