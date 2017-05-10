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
    
    
/** @file Bspline3DSetTemp.h
 * @brief Define template Bspline3DSet<bool ORTHO, bool TRUNC>
 *
 * Use specialization
 * Bspline3DSet<true,false> : orthorhombic unit cell
 * Bspline3DSet<false,false> : non-orthorhombic unit cell
 * Bspline3DSet<true,true> : orthorhombic unit cell with localized orbitals
 * Bspline3DSet<false,true> : non-orthorhombic unit cell with localized orbitals
 */
#ifndef QMCPLUSPLUS_BSPLINE3D_TEMPLATED_H
#define QMCPLUSPLUS_BSPLINE3D_TEMPLATED_H

#include "QMCWaveFunctions/Bspline3DSetBase.h"

namespace qmcplusplus
{

template<bool ORTHO, bool TRUNC>
class Bspline3DSet: public Bspline3DSetBase {};

/** specialized for Orthorhombic cell and no truncation*/
template<>
class Bspline3DSet<true,false>: public Bspline3DSetBase
{
public:

  Bspline3DSet() { }
  ~Bspline3DSet() { }

  inline void evaluate(const ParticleSet& e, int iat, ValueVector_t& vals)
  {
    bKnots.Find(e.R[iat][0],e.R[iat][1],e.R[iat][2]);
    for(int j=0; j<OrbitalSetSize; j++)
      vals[j]=bKnots.evaluate(*P[j]);
  }

  inline void evaluate(const ParticleSet& e, int iat,
                       ValueVector_t& vals, GradVector_t& grads, ValueVector_t& laps)
  {
    bKnots.FindAll(e.R[iat][0],e.R[iat][1],e.R[iat][2]);
    for(int j=0; j<OrbitalSetSize; j++)
      vals[j]=bKnots.evaluate(*P[j],grads[j],laps[j]);
  }

  inline void evaluate_notranspose(const ParticleSet& e, int first, int last,
                                   ValueMatrix_t& vals, GradMatrix_t& grads, ValueMatrix_t& laps)
  {
    for(int iat=first,i=0; iat<last; iat++,i++)
    {
      bKnots.FindAll(e.R[iat][0],e.R[iat][1],e.R[iat][2]);
      for(int j=0; j<OrbitalSetSize; j++)
        vals(i,j)=bKnots.evaluate(*P[j],grads(i,j),laps(i,j));
    }
  }
  void evaluate_notranspose(const ParticleSet& P, int first, int last
                            , ValueMatrix_t& logdet, GradMatrix_t& dlogdet, HessMatrix_t& grad_grad_logdet)
  {
    APP_ABORT("Need specialization of evaluate_notranspose() for grad_grad_logdet. \n");
  }
};

/** specialized for non-Orthorhombic cell and no truncation
*/
template<>
class Bspline3DSet<false,false>: public Bspline3DSetBase
{
public:

  Bspline3DSet()
  {
    Orthorhombic=false;
  }
  ~Bspline3DSet() { }

  inline void evaluate(const ParticleSet& e, int iat, ValueVector_t& vals)
  {
    PosType ru(Lattice.toUnit(e.R[iat]));
    bKnots.Find(ru[0],ru[1],ru[2]);
    for(int j=0; j<OrbitalSetSize; j++)
      vals[j]=bKnots.evaluate(*P[j]);
  }

  inline void evaluate(const ParticleSet& e, int iat,
                       ValueVector_t& vals, GradVector_t& grads, ValueVector_t& laps)
  {
    PosType ru(Lattice.toUnit(e.R[iat]));
    TinyVector<ValueType,3> gu;
    Tensor<ValueType,3> hess;
    bKnots.FindAll(ru[0],ru[1],ru[2]);
    for(int j=0; j<OrbitalSetSize; j++)
    {
      vals[j]=bKnots.evaluate(*P[j],gu,hess);
      grads[j]=dot(Lattice.G,gu);
      laps[j]=trace(hess,GGt);
    }
  }

  inline void evaluate_notranspose(const ParticleSet& e, int first, int last,
                                   ValueMatrix_t& vals, GradMatrix_t& grads, ValueMatrix_t& laps)
  {
    for(int iat=first,i=0; iat<last; iat++,i++)
    {
      PosType ru(Lattice.toUnit(e.R[iat]));
      TinyVector<ValueType,3> gu;
      Tensor<ValueType,3> hess;
      bKnots.FindAll(ru[0],ru[1],ru[2]);
      for(int j=0; j<OrbitalSetSize; j++)
      {
        vals(i,j)=bKnots.evaluate(*P[j],gu,hess);
        grads(i,j)=dot(Lattice.G,gu);
        laps(i,j)=trace(hess,GGt);
      }
    }
  }

  void evaluate_notranspose(const ParticleSet& P, int first, int last
                            , ValueMatrix_t& logdet, GradMatrix_t& dlogdet, HessMatrix_t& grad_grad_logdet)
  {
    APP_ABORT("Need specialization of evaluate_notranspose() for grad_grad_logdet. \n");
  }

};

/** specialized for Orthorhombic cell and no truncation*/
template<>
class Bspline3DSet<true,true>: public Bspline3DSetBase
{
public:
  Bspline3DSet()
  {
    Rcut2=1e6;
  }
  ~Bspline3DSet() { }

  inline void evaluate(const ParticleSet& e, int iat, ValueVector_t& vals)
  {
    PosType r(e.R[iat]);
    bKnots.Find(r[0],r[1],r[2]);
    for(int j=0; j<Centers.size(); j++)
    {
      if(bKnots.getSep2(r[0]-Centers[j][0],r[1]-Centers[j][1],r[2]-Centers[j][2])>Rcut2)
        vals[j]=0.0;//numeric_limits<T>::epsilon();
      else
        vals[j]=bKnots.evaluate(*P[j]);
    }
  }

  inline void evaluate(const ParticleSet& e, int iat,
                       ValueVector_t& vals, GradVector_t& grads, ValueVector_t& laps)
  {
    PosType r(e.R[iat]);
    bKnots.FindAll(r[0],r[1],r[2]);
    for(int j=0; j<Centers.size(); j++)
    {
      if(bKnots.getSep2(r[0]-Centers[j][0],r[1]-Centers[j][1],r[2]-Centers[j][2])>Rcut2)
      {
        vals[j]=0.0;//numeric_limits<T>::epsilon();
        grads[j]=0.0;
        laps[j]=0.0;
      }
      else
        vals[j]=bKnots.evaluate(*P[j],grads[j],laps[j]);
    }
  }

  inline void evaluate(const ParticleSet& e, int first, int last,
                       ValueMatrix_t& vals, GradMatrix_t& grads, ValueMatrix_t& laps)
  {
    for(int iat=first,i=0; iat<last; iat++,i++)
    {
      PosType r(e.R[iat]);
      bKnots.FindAll(r[0],r[1],r[2]);
      for(int j=0; j<Centers.size(); j++)
      {
        if(bKnots.getSep2(r[0]-Centers[j][0],r[1]-Centers[j][1],r[2]-Centers[j][2])>Rcut2)
        {
          vals(j,i)=0.0; //numeric_limits<T>::epsilon();
          grads(i,j)=0.0;
          laps(i,j)=0.0;
        }
        else
          vals(j,i)=bKnots.evaluate(*P[j],grads(i,j),laps(i,j));
      }
    }
  }

  void evaluate_notranspose(const ParticleSet& P, int first, int last
                            , ValueMatrix_t& logdet, GradMatrix_t& dlogdet, ValueMatrix_t& d2logdet)
  {
    APP_ABORT("Need specialization of evaluate_notranspose() for grad_grad_logdet. \n");
  }


  void evaluate_notranspose(const ParticleSet& P, int first, int last
                            , ValueMatrix_t& logdet, GradMatrix_t& dlogdet, HessMatrix_t& grad_grad_logdet)
  {
    APP_ABORT("Need specialization of evaluate_notranspose() for grad_grad_logdet. \n");
  }

};

/** specialized for non-Orthorhombic cell and truncation*/
template<>
class Bspline3DSet<false,true>: public Bspline3DSetBase
{
public:
  /** default constructure
   *
   * Set Rcut2 to a large number so that everything counts
   */
  Bspline3DSet()
  {
    Orthorhombic=false;
    Rcut2=1e6;
  }

  ~Bspline3DSet() { }

  inline void evaluate(const ParticleSet& e, int iat, ValueVector_t& vals)
  {
    PosType r(e.R[iat]);
    PosType ru(Lattice.toUnit(r));
    bKnots.Find(ru[0],ru[1],ru[2]);
    for(int j=0; j<Centers.size(); j++)
    {
      if(bKnots.getSep2(r[0]-Centers[j][0],r[1]-Centers[j][1],r[2]-Centers[j][2])>Rcut2)
        vals[j]=0.0;//numeric_limits<T>::epsilon();
      else
        vals[j]=bKnots.evaluate(*P[j]);
    }
  }

  inline void evaluate(const ParticleSet& e, int iat,
                       ValueVector_t& vals, GradVector_t& grads, ValueVector_t& laps)
  {
    PosType r(e.R[iat]);
    PosType ru(Lattice.toUnit(r));
    bKnots.FindAll(ru[0],ru[1],ru[2]);
    TinyVector<ValueType,3> gu;
    Tensor<ValueType,3> hess;
    for(int j=0; j<Centers.size(); j++)
    {
      if(bKnots.getSep2(r[0]-Centers[j][0],r[1]-Centers[j][1],r[2]-Centers[j][2])>Rcut2)
      {
        vals[j]=0.0;//numeric_limits<T>::epsilon();
        grads[j]=0.0;
        laps[j]=0.0;
      }
      else
      {
        vals[j]=bKnots.evaluate(*P[j],gu,hess);
        grads[j]=dot(Lattice.G,gu);
        laps[j]=trace(hess,GGt);
      }
    }
  }

  inline void evaluate_notranspose(const ParticleSet& e, int first, int last,
                                   ValueMatrix_t& vals, GradMatrix_t& grads, ValueMatrix_t& laps)
  {
    for(int iat=first,i=0; iat<last; iat++,i++)
    {
      PosType r(e.R[iat]);
      PosType ru(Lattice.toUnit(r));
      bKnots.FindAll(ru[0],ru[1],ru[2]);
      TinyVector<ValueType,3> gu;
      Tensor<ValueType,3> hess;
      for(int j=0; j<Centers.size(); j++)
      {
        if(bKnots.getSep2(r[0]-Centers[j][0],r[1]-Centers[j][1],r[2]-Centers[j][2])>Rcut2)
        {
          vals(i,j)=0.0; //numeric_limits<T>::epsilon();
          grads(i,j)=0.0;
          laps(i,j)=0.0;
        }
        else
        {
          vals(i,j)=bKnots.evaluate(*P[j],gu,hess);
          grads(i,j)=dot(Lattice.G,gu);
          laps(i,j)=trace(hess,GGt);
          //vals(j,i)=bKnots.evaluate(*P[j],grads(i,j),laps(i,j));
        }
      }
    }
  }

  void evaluate_notranspose(const ParticleSet& P, int first, int last
                            , ValueMatrix_t& logdet, GradMatrix_t& dlogdet, HessMatrix_t& grad_grad_logdet)
  {
    APP_ABORT("Need specialization of evaluate_notranspose() for grad_grad_logdet. \n");
  }

};
}
#endif
