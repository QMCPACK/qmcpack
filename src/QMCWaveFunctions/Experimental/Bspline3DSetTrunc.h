//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Miguel Morales, moralessilva2@llnl.gov, Lawrence Livermore National Laboratory
//                    Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//
// File created by: Miguel Morales, moralessilva2@llnl.gov, Lawrence Livermore National Laboratory
//////////////////////////////////////////////////////////////////////////////////////
    
    
/** @file Bspline3DSet.h
 * @brief Define Bspline3DSetBase and its derived classes
 *
 * - Bspline3DSet_Ortho : orthorhombic unit cell
 * - Bspline3DSet_Gen : non-orthorhombic unit cell
 * - Bspline3DSet_Ortho_Trunc: orthorhombic unit cell with localized orbitals
 * - Bspline3DSet_Gen_Trunc : non-orthorhombic unit cell with localized orbitals
 */
#ifndef QMCPLUSPLUS_BSPLINE3DSET_TRUCATEDANDTRANSLATED_H
#define QMCPLUSPLUS_BSPLINE3DSET_TRUCATEDANDTRANSLATED_H

#include "QMCWaveFunctions/Bspline3DSetBase.h"

namespace qmcplusplus
{

/** Specialized for Maximally Localized Wavetions
 *
 * Grids for localized orbitals are truncated.
 */
struct Bspline3DSet_MLW: public Bspline3DSetBase
{
  Bspline3DSet_MLW();
  ~Bspline3DSet_MLW();

  RealType Lx,LxInv;
  RealType Ly,LyInv;
  RealType Lz,LzInv;
  RealType LxSq, LySq, LzSq;

  inline PosType translate(const PosType& r, int j)
  {
    PosType rtr=r-Origins[j]; // first shift it by the origin
    rtr[0]-=std::floor(rtr[0]*Lattice.G[0])*Lattice.R[0];
    rtr[1]-=std::floor(rtr[1]*Lattice.G[4])*Lattice.R[4];
    rtr[2]-=std::floor(rtr[2]*Lattice.G[8])*Lattice.R[8];
    //rtr[0]-=std::floor(rtr[0]*LxInv)*Lx;
    //rtr[1]-=std::floor(rtr[1]*LyInv)*Ly;
    //rtr[2]-=std::floor(rtr[2]*LzInv)*Lz;
    return rtr;
  }

  /* return the distance between the center with PBC */
  inline RealType getSep2(RealType x, RealType y, RealType z)
  {
    //x-=nearbyint(x*Lattice.G[0])*Lattice.R[0];
    //y-=nearbyint(y*Lattice.G[4])*Lattice.R[4];
    //z-=nearbyint(z*Lattice.G[8])*Lattice.R[8];
    RealType x1=std::fmod(x*Lattice.G[0],1.0);
    x=Lattice.R[0]*(x1-static_cast<int>(2.0*x1));
    RealType y1=std::fmod(y*Lattice.G[4],1.0);
    y=Lattice.R[4]*(y1-static_cast<int>(2.0*y1));
    RealType z1=std::fmod(z*Lattice.G[8],1.0);
    z=Lattice.R[8]*(z1-static_cast<int>(2.0*z1));
    return x*x+y*y+z*z;
  }

  void evaluate(const ParticleSet& P, int iat, ValueVector_t& psi);
  void evaluate(const ParticleSet& P, int iat,
                ValueVector_t& psi, GradVector_t& dpsi, ValueVector_t& d2psi);
  void evaluate(const ParticleSet& P, int iat,
                ValueVector_t& psi, GradVector_t& dpsi, HessVector_t& gg_psi)
  {
    APP_ABORT("Need specialization of evaluate(iat) for HessVector. \n");
  }
  void evaluate(const ParticleSet& P, int first, int last,
                ValueMatrix_t& logdet, GradMatrix_t& dlogdet, ValueMatrix_t& d2logdet);
  void evaluate_notranspose(const ParticleSet& P, int first, int last,
                            ValueMatrix_t& logdet, GradMatrix_t& dlogdet, ValueMatrix_t& d2logdet);
  void evaluate_notranspose(const ParticleSet& P, int first, int last
                            , ValueMatrix_t& logdet, GradMatrix_t& dlogdet, HessMatrix_t& grad_grad_logdet)
  {
    APP_ABORT("Need specialization of evaluate_notranspose() for grad_grad_logdet. \n");
  }

};

//  /** Specialized for non-Orthorhombic cell no truncation*/
//  struct Bspline3DSet_Gen_Trunc: public Bspline3DSetBase
//  {
//
//    Bspline3DSet_Gen_Trunc() {Orthorhombic=false;}
//    ~Bspline3DSet_Gen_Trunc() { }
//
//    void evaluate(const ParticleSet& P, int iat, ValueVector_t& psi);
//    void evaluate(const ParticleSet& P, int iat,
//        ValueVector_t& psi, GradVector_t& dpsi, ValueVector_t& d2psi);
//    void evaluate(const ParticleSet& P, int first, int last,
//        ValueMatrix_t& logdet, GradMatrix_t& dlogdet, ValueMatrix_t& d2logdet);
//  };

}
#endif
