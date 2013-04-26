/////////////////////////////////////////////////////////////////
// (c) Copyright 2007-  Jeongnim Kim
//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////
//   Modified by Jeongnim Kim for qmcpack
//   National Center for Supercomputing Applications &
//   Materials Computation Center
//   University of Illinois, Urbana-Champaign
//   Urbana, IL 61801
//   e-mail: jnkim@ncsa.uiuc.edu
//
// Supported by
//   National Center for Supercomputing Applications, UIUC
//   Materials Computation Center, UIUC
//////////////////////////////////////////////////////////////////
// -*- C++ -*-
/** @file Bspline3DSet.h
 * @brief Define Bspline3DSetBase and its derived classes
 *
 * - Bspline3DSet_Ortho : orthorhombic unit cell
 * - Bspline3DSet_Gen : non-orthorhombic unit cell
 * - Bspline3DSet_Ortho_Trunc: orthorhombic unit cell with localized orbitals
 * - Bspline3DSet_Gen_Trunc : non-orthorhombic unit cell with localized orbitals
 */
#ifndef QMCPLUSPLUS_BSPLINE3DSET_H
#define QMCPLUSPLUS_BSPLINE3DSET_H

#include "QMCWaveFunctions/Bspline3DSetBase.h"

namespace qmcplusplus
{

/** Specialized for Orthorhombic cell and no truncation*/
struct Bspline3DSet_Ortho: public Bspline3DSetBase
{
  Bspline3DSet_Ortho() {}
  ~Bspline3DSet_Ortho() { }

  void evaluate(const ParticleSet& P, int iat, ValueVector_t& psi);
  void evaluate(const ParticleSet& P, int iat,
                ValueVector_t& psi, GradVector_t& dpsi, ValueVector_t& d2psi);
  void evaluate(const ParticleSet& P, int iat,
                ValueVector_t& psi, GradVector_t& dpsi, HessVector_t& gg_psi)
  {
    APP_ABORT("Need specialization of evaluate(iat) for HessVector. \n");
  }

  void evaluate_notranspose(const ParticleSet& P, int first, int last,
                            ValueMatrix_t& logdet, GradMatrix_t& dlogdet, ValueMatrix_t& d2logdet);
  void evaluate_notranspose(const ParticleSet& P, int first, int last
                            , ValueMatrix_t& logdet, GradMatrix_t& dlogdet, HessMatrix_t& grad_grad_logdet)
  {
    APP_ABORT("Need specialization of evaluate_notranspose() for grad_grad_logdet. \n");
  }

  SPOSetBase* makeClone() const;

};

/** Specialized for Orthorhombic cell and no truncation*/
struct Bspline3DSet_Gen: public Bspline3DSetBase
{

  Bspline3DSet_Gen()
  {
    Orthorhombic=false;
  }
  ~Bspline3DSet_Gen() { }

  void evaluate(const ParticleSet& P, int iat, ValueVector_t& psi);
  void evaluate(const ParticleSet& P, int iat,
                ValueVector_t& psi, GradVector_t& dpsi, ValueVector_t& d2psi);
  void evaluate(const ParticleSet& P, int iat,
                ValueVector_t& psi, GradVector_t& dpsi, HessVector_t& gg_psi)
  {
    APP_ABORT("Need specialization of evaluate(iat) for HessVector. \n");
  }
  void evaluate_notranspose(const ParticleSet& P, int first, int last,
                            ValueMatrix_t& logdet, GradMatrix_t& dlogdet, ValueMatrix_t& d2logdet);
  void evaluate_notranspose(const ParticleSet& P, int first, int last
                            , ValueMatrix_t& logdet, GradMatrix_t& dlogdet, HessMatrix_t& grad_grad_logdet)
  {
    APP_ABORT("Need specialization of evaluate_notranspose() for grad_grad_logdet. \n");
  }

  SPOSetBase* makeClone() const;
};

/** Specialized for Orthorhombic cell and truncation*/
struct Bspline3DSet_Ortho_Trunc: public Bspline3DSetBase
{
  Bspline3DSet_Ortho_Trunc() {}
  ~Bspline3DSet_Ortho_Trunc() { }

  void evaluate(const ParticleSet& P, int iat, ValueVector_t& psi);
  void evaluate(const ParticleSet& P, int iat,
                ValueVector_t& psi, GradVector_t& dpsi, ValueVector_t& d2psi);
  void evaluate(const ParticleSet& P, int iat,
                ValueVector_t& psi, GradVector_t& dpsi, HessVector_t& gg_psi)
  {
    APP_ABORT("Need specialization of evaluate(iat) for HessVector. \n");
  }
  void evaluate_notranspose(const ParticleSet& P, int first, int last,
                            ValueMatrix_t& logdet, GradMatrix_t& dlogdet, ValueMatrix_t& d2logdet);
  void evaluate_notranspose(const ParticleSet& P, int first, int last
                            , ValueMatrix_t& logdet, GradMatrix_t& dlogdet, HessMatrix_t& grad_grad_logdet)
  {
    APP_ABORT("Need specialization of evaluate_notranspose() for grad_grad_logdet. \n");
  }

  SPOSetBase* makeClone() const;
};

/** Specialized for non-Orthorhombic cell no truncation*/
struct Bspline3DSet_Gen_Trunc: public Bspline3DSetBase
{

  Bspline3DSet_Gen_Trunc()
  {
    Orthorhombic=false;
  }
  ~Bspline3DSet_Gen_Trunc() { }

  void evaluate(const ParticleSet& P, int iat, ValueVector_t& psi);
  void evaluate(const ParticleSet& P, int iat,
                ValueVector_t& psi, GradVector_t& dpsi, ValueVector_t& d2psi);
  void evaluate(const ParticleSet& P, int iat,
                ValueVector_t& psi, GradVector_t& dpsi, HessVector_t& gg_psi)
  {
    APP_ABORT("Need specialization of evaluate(iat) for HessVector. \n");
  }
  void evaluate_notranspose(const ParticleSet& P, int first, int last,
                            ValueMatrix_t& logdet, GradMatrix_t& dlogdet, ValueMatrix_t& d2logdet);
  void evaluate_notranspose(const ParticleSet& P, int first, int last
                            , ValueMatrix_t& logdet, GradMatrix_t& dlogdet, HessMatrix_t& grad_grad_logdet)
  {
    APP_ABORT("Need specialization of evaluate_notranspose() for grad_grad_logdet. \n");
  }

  SPOSetBase* makeClone() const;
};

#if defined(QMC_COMPLEX)
struct Bspline3DSet_Twist: public Bspline3DSetBase
{
  Bspline3DSet_Twist()
  {
    Orthorhombic=false;
  }
  ~Bspline3DSet_Twist() { }

  void evaluate(const ParticleSet& P, int iat, ValueVector_t& psi);
  void evaluate(const ParticleSet& P, int iat,
                ValueVector_t& psi, GradVector_t& dpsi, ValueVector_t& d2psi);
  void evaluate(const ParticleSet& P, int iat,
                ValueVector_t& psi, GradVector_t& dpsi, HessVector_t& gg_psi)
  {
    APP_ABORT("Need specialization of evaluate(iat) for HessVector. \n");
  }
  void evaluate_notranspose(const ParticleSet& P, int first, int last,
                            ValueMatrix_t& logdet, GradMatrix_t& dlogdet, ValueMatrix_t& d2logdet);
  void evaluate_notranspose(const ParticleSet& P, int first, int last
                            , ValueMatrix_t& logdet, GradMatrix_t& dlogdet, HessMatrix_t& grad_grad_logdet)
  {
    APP_ABORT("Need specialization of evaluate_notranspose() for grad_grad_logdet. \n");
  }

  SPOSetBase* makeClone() const;
};
#endif

}
#endif
/***************************************************************************
 * $RCSfile$   $Author: jnkim $
 * $Revision: 2013 $   $Date: 2007-05-22 16:47:09 -0500 (Tue, 22 May 2007) $
 * $Id: TricubicBsplineSet.h 2013 2007-05-22 21:47:09Z jnkim $
 ***************************************************************************/
