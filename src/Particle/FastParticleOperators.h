//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////
    
    



/** @file FastParticleOperators.h
 * @brief template functions to support conversion and handling of boundary conditions.
 */
#ifndef OHMMS_FAST_PARTICLE_OPERATORS_H
#define OHMMS_FAST_PARTICLE_OPERATORS_H

#include <simd/simd.hpp>

namespace qmcplusplus
{
/** Dummy template class to be specialized
 *
 * - T1 the datatype to be transformed
 * - T2 the transformation matrix
 * - ORTHO true, if only Diagonal Elements are used
 */
template<class T1, class T2, unsigned D, bool ORTHO> struct ConvertPosUnit { };

/** Specialized ConvertPosUnit for ParticleAttrib<T,3>, Tensor<T,3> and true
 */
template<class T>
struct ConvertPosUnit<ParticleAttrib<TinyVector<T,3> >,Tensor<T,3>, 3, true>
{

  typedef ParticleAttrib<TinyVector<T,3> > Array_t;
  typedef Tensor<T,3>                      Transformer_t;

  /** apply the transformation matrix, pout[i] = dot(pin[i],X) for i=[first,last)
   *
   * @param pin input data Array to be transformed
   * @param X  transformation matrix which operates to the left vector
   * @param pout outout data Array
   * @param first the first index
   * @param last the last index
   */
  inline static void
  apply(const Array_t& pin, const Transformer_t& X, Array_t& pout, int first, int last)
  {
    register T xx=X[0], yy=X[4], zz=X[8];
#pragma ivdep
    for(int i=first; i<last; i++)
    {
      pout[i][0]=pin[i][0]*xx;
      pout[i][1]=pin[i][1]*yy;
      pout[i][2]=pin[i][2]*zz;
    }
  }

  /** apply the transformation matrix, pout[i] = dot(X,pin[i]) for i=[first,last)
   *
   * @param X  transformation matrix which operates to the right vector
   * @param pin input data Array to be transformed
   * @param pout outout data Array
   * @param first the first index
   * @param last the last index
   */
  inline static void
  apply(const Transformer_t& X, const Array_t& pin,  Array_t& pout, int first, int last)
  {
    register T xx=X[0], yy=X[4], zz=X[8];
#pragma ivdep
    for(int i=first; i<last; i++)
    {
      pout[i][0]=pin[i][0]*xx;
      pout[i][1]=pin[i][1]*yy;
      pout[i][2]=pin[i][2]*zz;
    }
  }

  /** apply the transformation matrix, pinout[i] = dot(pinout[i],X) for i=[first,last)
   *
   * @param pinout input/output data Array to be transformed
   * @param X  transformation matrix which operates to the left vector
   * @param first the first index
   * @param last the last index
   */
  inline static void
  apply(Array_t& pinout, const Transformer_t& X,int first, int last)
  {
    register T xx=X[0], yy=X[4], zz=X[8];
#pragma ivdep
    for(int i=first; i<last; i++)
    {
      pinout[i][0]*=xx;
      pinout[i][1]*=yy;
      pinout[i][2]*=zz;
    }
  }

  /** apply the transformation matrix, pinout[i] = dot(X,pinout[i]) for i=[first,last)
   *
   * @param X  transformation matrix which operates to the right vector
   * @param pinout input/output data Array to be transformed
   * @param first the first index
   * @param last the last index
   */
  inline static void
  apply(const Transformer_t& X, Array_t& pinout, int first, int last)
  {
    register T xx=X[0], yy=X[4], zz=X[8];
#pragma ivdep
    for(int i=first; i<last; i++)
    {
      pinout[i][0]*=xx;
      pinout[i][1]*=yy;
      pinout[i][2]*=zz;
    }
  }
};

template<class T>
struct ConvertPosUnit<ParticleAttrib<TinyVector<T,3> >,Tensor<T,3>, 3, false>
{

  typedef ParticleAttrib<TinyVector<T,3> > Array_t;
  typedef Tensor<T,3>                      Transformer_t;

  inline static void
  apply(const Array_t& pin, const Transformer_t& X, Array_t& pout, int first, int last)
  {
    register T x00=X[0],x01=X[1],x02=X[2],
               x10=X[3],x11=X[4],x12=X[5],
               x20=X[6],x21=X[7],x22=X[8];
#pragma ivdep
    for(int i=first; i<last; i++)
    {
      pout[i][0]=pin[i][0]*x00+pin[i][1]*x10+pin[i][2]*x20;
      pout[i][1]=pin[i][0]*x01+pin[i][1]*x11+pin[i][2]*x21;
      pout[i][2]=pin[i][0]*x02+pin[i][1]*x12+pin[i][2]*x22;
    }
  }

  inline static void
  apply(const Transformer_t& X, const Array_t& pin,  Array_t& pout, int first, int last)
  {
    register T x00=X[0],x01=X[1],x02=X[2],
               x10=X[3],x11=X[4],x12=X[5],
               x20=X[6],x21=X[7],x22=X[8];
#pragma ivdep
    for(int i=first; i<last; i++)
    {
      pout[i][0]=pin[i][0]*x00+pin[i][1]*x01+pin[i][2]*x02;
      pout[i][1]=pin[i][0]*x10+pin[i][1]*x11+pin[i][2]*x12;
      pout[i][2]=pin[i][0]*x20+pin[i][1]*x21+pin[i][2]*x22;
    }
  }

  inline static void
  apply(Array_t& pinout, const Transformer_t& X,int first, int last)
  {
    register T x00=X[0],x01=X[1],x02=X[2],
               x10=X[3],x11=X[4],x12=X[5],
               x20=X[6],x21=X[7],x22=X[8];
#pragma ivdep
    for(int i=first; i<last; i++)
    {
      T _x(pinout[i][0]),_y(pinout[i][1]),_z(pinout[i][2]);
      pinout[i][0]=_x*x00+_y*x10+_z*x20;
      pinout[i][1]=_x*x01+_y*x11+_z*x21;
      pinout[i][2]=_x*x02+_y*x12+_z*x22;
    }
  }

  inline static void
  apply(const Transformer_t& X, Array_t& pinout, int first, int last)
  {
    register T x00=X[0],x01=X[1],x02=X[2],
               x10=X[3],x11=X[4],x12=X[5],
               x20=X[6],x21=X[7],x22=X[8];
#pragma ivdep
    for(int i=first; i<last; i++)
    {
      T _x(pinout[i][0]),_y(pinout[i][1]),_z(pinout[i][2]);
      pinout[i][0]=_x*x00+_y*x01+_z*x02;
      pinout[i][1]=_x*x10+_y*x11+_z*x12;
      pinout[i][2]=_x*x20+_y*x21+_z*x22;
    }
  }
};


/** Specialized ConvertPosUnit for ParticleAttrib<T,2>, Tensor<T,2> and true
 */
template<class T>
struct ConvertPosUnit<ParticleAttrib<TinyVector<T,2> >,Tensor<T,2>, 2, true>
{

  typedef ParticleAttrib<TinyVector<T,2> > Array_t;
  typedef Tensor<T,2>                      Transformer_t;

  inline static void
  apply(const Array_t& pin, const Transformer_t& X, Array_t& pout, int first, int last)
  {
    register T xx=X[0], yy=X[3];
#pragma ivdep
    for(int i=first; i<last; i++)
    {
      pout[i][0]=pin[i][0]*xx;
      pout[i][1]=pin[i][1]*yy;
    }
  }

  inline static void
  apply(const Transformer_t& X, const Array_t& pin,  Array_t& pout, int first, int last)
  {
    register T xx=X[0], yy=X[3];
#pragma ivdep
    for(int i=first; i<last; i++)
    {
      pout[i][0]=pin[i][0]*xx;
      pout[i][1]=pin[i][1]*yy;
    }
  }

  inline static void
  apply(Array_t& pinout, const Transformer_t& X,int first, int last)
  {
    register T xx=X[0], yy=X[3];
#pragma ivdep
    for(int i=first; i<last; i++)
    {
      pinout[i][0]*=xx;
      pinout[i][1]*=yy;
    }
  }

  inline static void
  apply(const Transformer_t& X, Array_t& pinout, int first, int last)
  {
    register T xx=X[0], yy=X[3];
#pragma ivdep
    for(int i=first; i<last; i++)
    {
      pinout[i][0]*=xx;
      pinout[i][1]*=yy;
    }
  }
};

template<class T>
struct ConvertPosUnit<ParticleAttrib<TinyVector<T,2> >,Tensor<T,2>, 2, false>
{

  typedef ParticleAttrib<TinyVector<T,2> > Array_t;
  typedef Tensor<T,2>                      Transformer_t;

  inline static void
  apply(const Array_t& pin, const Transformer_t& X, Array_t& pout, int first, int last)
  {
    register T x00=X[0],x01=X[1], x10=X[2],x11=X[3];
#pragma ivdep
    for(int i=first; i<last; i++)
    {
      pout[i][0]=pin[i][0]*x00+pin[i][1]*x10;
      pout[i][1]=pin[i][0]*x01+pin[i][1]*x11;
    }
  }

  inline static void
  apply(const Transformer_t& X, const Array_t& pin,  Array_t& pout, int first, int last)
  {
    register T x00=X[0],x01=X[1],x10=X[2],x11=X[3];
#pragma ivdep
    for(int i=first; i<last; i++)
    {
      pout[i][0]=pin[i][0]*x00+pin[i][1]*x01;
      pout[i][1]=pin[i][0]*x10+pin[i][1]*x11;
    }
  }

  inline static void
  apply(Array_t& pinout, const Transformer_t& X,int first, int last)
  {
    register T x00=X[0],x01=X[1],x10=X[2],x11=X[3];
#pragma ivdep
    for(int i=first; i<last; i++)
    {
      T _x(pinout[i][0]),_y(pinout[i][1]);
      pinout[i][0]=_x*x00+_y*x10;
      pinout[i][1]=_x*x01+_y*x11;
    }
  }

  inline static void
  apply(const Transformer_t& X, Array_t& pinout, int first, int last)
  {
    register T x00=X[0],x01=X[1],x10=X[2],x11=X[3];
#pragma ivdep
    for(int i=first; i<last; i++)
    {
      T _x(pinout[i][0]),_y(pinout[i][1]);
      pinout[i][0]=_x*x00+_y*x01;
      pinout[i][1]=_x*x10+_y*x11;
    }
  }
};

#define SUPERCELL_BOUNDARY_LIMITS(T)               \
  const T epsilon = -std::numeric_limits<T>::epsilon(); \
  const T plus_one = 1.0


#define THREE_DIM_BOUNDARY_BLOCK(X,Y,Z,EPS,PLUSONE) \
    if(X<EPS)           X+=PLUSONE;                 \
    else if(X>=PLUSONE) X-=PLUSONE;                 \
    if(Y<EPS)           Y+=PLUSONE;                 \
    else if(Y>=PLUSONE) Y-=PLUSONE;                 \
    if(Z<EPS)           Z+=PLUSONE;                 \
    else if(Z>=PLUSONE) Z-=PLUSONE

#define TWO_DIM_BOUNDARY_BLOCK(X,Y,EPS,PLUSONE) \
    if(X<EPS)           X+=PLUSONE;                 \
    else if(X>=PLUSONE) X-=PLUSONE;                 \
    if(Y<EPS)           Y+=PLUSONE;                 \
    else if(Y>=PLUSONE) Y-=PLUSONE


/** Dummy template class to apply boundary conditions */
template<class T1, class T2, unsigned D, bool ORTHO> struct ApplyBConds { };

template<class T>
struct ApplyBConds<ParticleAttrib<TinyVector<T,3> >, Tensor<T,3>, 3,true>
{

  typedef ParticleAttrib<TinyVector<T,3> > Array_t;
  typedef typename Array_t::Type_t         Component_t;
  typedef Tensor<T,3>                      Transformer_t;

  /** Apply boundary condition on the LatticeUnit vectors and return LatticeUnit vectors
   * @param pin input array in the LatticeUnit Unit
   * @param pout output array
   * @param first starting index
   * @param last ending index
   *
   * Move the components to [0,1)
   */
  inline static void
  Unit2Unit(const Array_t& pin, Array_t& pout, int first, int last)
  {
    SUPERCELL_BOUNDARY_LIMITS(T);
#pragma ivdep
    for(int i=first; i<last; i++)
    {
      T x(pin[i][0]),y(pin[i][1]),z(pin[i][2]);
      THREE_DIM_BOUNDARY_BLOCK(x,y,z,epsilon,plus_one);
      pout[i][0]=x;
      pout[i][1]=y;
      pout[i][2]=z;
    }
  }

  /** Apply boundary condition on the LatticeUnit vectors and return Cartesian vectors
   * @param pin input array in the LatticeUnit
   * @param R transformation matrix from LatticeUnit to CartesianUnit
   * @param pout input array in the CartesianUnit
   * @param first starting index
   * @param last ending index
   *
   * pout = dot(applybconds(pin),R)
   * - applybconds move all to [0,1)
   * - dot(X,R) convert to CartesianUnit
   */
  inline static void
  Unit2Cart(const Array_t& pin, const Transformer_t& R, Array_t& pout, int first, int last)
  {
    SUPERCELL_BOUNDARY_LIMITS(T);
    register T r00=R[0], r11=R[4], r22=R[8];
#pragma ivdep
    for(int i=first; i<last; i++)
    {
      T x(pin[i][0]),y(pin[i][1]),z(pin[i][2]);
      THREE_DIM_BOUNDARY_BLOCK(x,y,z,epsilon,plus_one);
      pout[i][0]=r00*x;
      pout[i][1]=r11*y;
      pout[i][2]=r22*z;
    }
  }


  /** Apply boundary condition on the Cartesin vectors and return LatticeUnit vectors
   * @param pin input array in the Cartesian Unit
   * @param G transformation matrix from CartesianUnit to LatticeUnit
   * @param pout input array in the LatticeUnit
   * @param first starting index
   * @param last ending index
   *
   * pout = applybconds(dot(pin,G))
   * - dot(pin,G) convert to LatticeUnit
   * - applybconds move all to [0,1)
   */
  inline static void
  Cart2Unit(const Array_t& pin, const Transformer_t& G, Array_t& pout, int first, int last)
  {
    SUPERCELL_BOUNDARY_LIMITS(T);
    register T g00=G[0], g11=G[4], g22=G[8];
#pragma ivdep
    for(int i=first; i<last; i++)
    {
      T x(pin[i][0]*g00),y(pin[i][1]*g11),z(pin[i][2]*g22);
      THREE_DIM_BOUNDARY_BLOCK(x,y,z,epsilon,plus_one);
      pout[i][0]=x;
      pout[i][1]=y;
      pout[i][2]=z;
    }
  }

  /** Apply boundary condition on the Cartesin vectors and return Cartesian vectors
   * @param pin input array in the Cartesian Unit
   * @param G transformation matrix from CartesianUnit to LatticeUnit
   * @param R transformation matrix from LatticeUnit to CartesianUnit
   * @param pout input array in the Cartesian Unit
   * @param first starting index
   * @param last ending index
   *
   * pout = dot(applybconds(dot(pin,G)),R)
   * - dot(pin,G) convert to LatticeUnit
   * - applybconds move all to [0,1)
   * - dot(X,R) convert to CartesianUnit
   */
  inline static void
  Cart2Cart(const Array_t& pin, const Transformer_t& G, const Transformer_t& R, Array_t& pout,
            int first, int last)
  {
    SUPERCELL_BOUNDARY_LIMITS(T);
    register T g00=G[0], g11=G[4], g22=G[8];
    register T r00=R[0], r11=R[4], r22=R[8];
#pragma ivdep
    for(int i=first; i<last; i++)
    {
      T x(pin[i][0]*g00),y(pin[i][1]*g11),z(pin[i][2]*g22);
      THREE_DIM_BOUNDARY_BLOCK(x,y,z,epsilon,plus_one);
      pout[i][0]=r00*x;
      pout[i][1]=r11*y;
      pout[i][2]=r22*z;
    }
  }

  /** Apply boundary condition on the LatticeUnit vectors and return LatticeUnit vectors
   * @param pinout input/output array in the LatticeUnit
   * @param first starting index
   * @param last ending index
   *
   * Move the components to [0,1)
   */
  inline static void
  Unit2Unit(Array_t& pinout, int first, int last)
  {
    SUPERCELL_BOUNDARY_LIMITS(T);
#pragma ivdep
    for(int i=first; i<last; i++)
    {
      T x(pinout[i][0]),y(pinout[i][1]),z(pinout[i][2]);
      THREE_DIM_BOUNDARY_BLOCK(x,y,z,epsilon,plus_one);
      pinout[i][0]=x;
      pinout[i][1]=y;
      pinout[i][2]=z;
    }
  }

  /** Apply boundary condition on the Cartesian vectors and return Cartesian vectors
   * @param pinout input/output array in the LatticeUnit
   * @param G transformation matrix from CartesianUnit to LatticeUnit
   * @param R transformation matrix from LatticeUnit to CartesianUnit
   * @param first starting index
   * @param last ending index
   *
   * pinout <- dot(applybconds(dot(pinout,G)),R)
   * - dot(pin,G) convert to LatticeUnit
   * - applybconds move all to [0,1)
   * - dot(X,R) convert to CartesianUnit
   */
  inline static void
  Cart2Cart(Array_t& pinout, const Transformer_t& G, const Transformer_t& R, int first, int last)
  {
    SUPERCELL_BOUNDARY_LIMITS(T);
    register T g00=G[0], g11=G[4], g22=G[8];
    register T r00=R[0], r11=R[4], r22=R[8];
#pragma ivdep
    for(int i=first; i<last; i++)
    {
      T x(pinout[i][0]*g00),y(pinout[i][1]*g11),z(pinout[i][2]*g22);
      THREE_DIM_BOUNDARY_BLOCK(x,y,z,epsilon,plus_one);
      pinout[i][0]=r00*x;
      pinout[i][1]=r11*y;
      pinout[i][2]=r22*z;
    }
  }

  static inline Component_t Unit2Unit(const Component_t& pos)
  {
    SUPERCELL_BOUNDARY_LIMITS(T);
    register T x(pos[0]),y(pos[1]),z(pos[2]);
    THREE_DIM_BOUNDARY_BLOCK(x,y,z,epsilon,plus_one);
    return Component_t(x,y,z);
  }

  static inline Component_t Cart2Unit(const Component_t& pos, const Transformer_t& G)
  {
    SUPERCELL_BOUNDARY_LIMITS(T);
    register T x(pos[0]*G[0]), y(pos[1]*G[4]),z(pos[2]*G[8]);
    THREE_DIM_BOUNDARY_BLOCK(x,y,z,epsilon,plus_one);
    return Component_t(x,y,z);
  }

  static inline Component_t Unit2Cart(const Component_t& pos, const Transformer_t& R)
  {
    SUPERCELL_BOUNDARY_LIMITS(T);
    register T x(pos[0]),y(pos[1]),z(pos[2]);
    THREE_DIM_BOUNDARY_BLOCK(x,y,z,epsilon,plus_one);
    return Component_t(x*R[0],y*R[4],z*R[8]);
  }

  static inline Component_t Cart2Cart(const Component_t& pos, const Transformer_t& G, const Transformer_t& R)
  {
    SUPERCELL_BOUNDARY_LIMITS(T);
    register T x(pos[0]*G[0]),y(pos[1]*G[4]),z(pos[2]*G[8]);
    THREE_DIM_BOUNDARY_BLOCK(x,y,z,epsilon,plus_one);
    return Component_t(x*R[0],y*R[4],z*R[8]);
  }
};

template<class T>
struct ApplyBConds<ParticleAttrib<TinyVector<T,3> >, Tensor<T,3>, 3,false>
{

  typedef ParticleAttrib<TinyVector<T,3> > Array_t;
  typedef typename Array_t::Type_t         Component_t;
  typedef Tensor<T,3>                      Transformer_t;

  inline static void
  Cart2Cart(const Array_t& pin, const Transformer_t& G, const Transformer_t& R, Array_t& pout, int first, int last)
  {
    SUPERCELL_BOUNDARY_LIMITS(T);
    register T g00=G[0], g01=G[1],g02=G[2], g10=G[3],g11=G[4],g12=G[5], g20=G[6],g21=G[7],g22=G[8];
    register T r00=R[0], r01=R[1],r02=R[2], r10=R[3],r11=R[4],r12=R[5], r20=R[6],r21=R[7],r22=R[8];
#pragma ivdep
    for(int i=first; i<last; i++)
    {
      T x(pin[i][0]*g00+pin[i][1]*g10+pin[i][2]*g20);
      T y(pin[i][0]*g01+pin[i][1]*g11+pin[i][2]*g21);
      T z(pin[i][0]*g02+pin[i][1]*g12+pin[i][2]*g22);
      THREE_DIM_BOUNDARY_BLOCK(x,y,z,epsilon,plus_one);
      pout[i][0]=x*r00+y*r10+z*r20;
      pout[i][1]=x*r01+y*r11+z*r21;
      pout[i][2]=x*r02+y*r12+z*r22;
    }
  }

  inline static void
  Cart2Unit(const Array_t& pin, const Transformer_t& G, Array_t& pout, int first, int last)
  {
    SUPERCELL_BOUNDARY_LIMITS(T);
    register T g00=G[0],g01=G[1],g02=G[2], g10=G[3],g11=G[4],g12=G[5], g20=G[6],g21=G[7],g22=G[8];
#pragma ivdep
    for(int i=first; i<last; i++)
    {
      T x(pin[i][0]*g00+pin[i][1]*g10+pin[i][2]*g20);
      T y(pin[i][0]*g01+pin[i][1]*g11+pin[i][2]*g21);
      T z(pin[i][0]*g02+pin[i][1]*g12+pin[i][2]*g22);
      THREE_DIM_BOUNDARY_BLOCK(x,y,z,epsilon,plus_one);
      pout[i][0]=x;
      pout[i][1]=y;
      pout[i][2]=z;
    }
  }

  inline static void
  Unit2Cart(const Array_t& pin, const Transformer_t& R, Array_t& pout, int first, int last)
  {
    SUPERCELL_BOUNDARY_LIMITS(T);
    register T r00=R[0],r01=R[1],r02=R[2], r10=R[3],r11=R[4],r12=R[5], r20=R[6],r21=R[7],r22=R[8];
#pragma ivdep
    for(int i=first; i<last; i++)
    {
      T x(pin[i][0]), y(pin[i][1]),z(pin[i][2]);
      THREE_DIM_BOUNDARY_BLOCK(x,y,z,epsilon,plus_one);
      pout[i][0]=x*r00+y*r10+z*r20;
      pout[i][1]=x*r01+y*r11+z*r21;
      pout[i][2]=x*r02+y*r12+z*r22;
    }
  }

  inline static void
  Unit2Unit(const Array_t& pin, Array_t& pout, int first, int last)
  {
    SUPERCELL_BOUNDARY_LIMITS(T);
#pragma ivdep
    for(int i=first; i<last; i++)
    {
      T x(pin[i][0]),y(pin[i][1]),z(pin[i][2]);
      THREE_DIM_BOUNDARY_BLOCK(x,y,z,epsilon,plus_one);
      pout[i][0]=x;
      pout[i][1]=y;
      pout[i][2]=z;
    }
  }

  inline static void
  Unit2Unit(Array_t& pinout, int first, int last)
  {
    SUPERCELL_BOUNDARY_LIMITS(T);
#pragma ivdep
    for(int i=first; i<last; i++)
    {
      T x(pinout[i][0]),y(pinout[i][1]),z(pinout[i][2]);
      THREE_DIM_BOUNDARY_BLOCK(x,y,z,epsilon,plus_one);
      pinout[i][0]=x;
      pinout[i][1]=y;
      pinout[i][2]=z;
    }
  }

  inline static void
  Cart2Cart(Array_t& pinout, const Transformer_t& G, const Transformer_t& R, int first, int last)
  {
    SUPERCELL_BOUNDARY_LIMITS(T);
    register T g00=G[0],g01=G[1],g02=G[2], g10=G[3],g11=G[4],g12=G[5], g20=G[6],g21=G[7],g22=G[8];
    register T r00=R[0],r01=R[1],r02=R[2], r10=R[3],r11=R[4],r12=R[5], r20=R[6],r21=R[7],r22=R[8];
#pragma ivdep
    for(int i=first; i<last; i++)
    {
      T x(pinout[i][0]*g00+pinout[i][1]*g10+pinout[i][2]*g20);
      T y(pinout[i][0]*g01+pinout[i][1]*g11+pinout[i][2]*g21);
      T z(pinout[i][0]*g02+pinout[i][1]*g12+pinout[i][2]*g22);
      THREE_DIM_BOUNDARY_BLOCK(x,y,z,epsilon,plus_one);
      pinout[i][0]=x*r00+y*r10+z*r20;
      pinout[i][1]=x*r01+y*r11+z*r21;
      pinout[i][2]=x*r02+y*r12+z*r22;
    }
  }

  static inline Component_t Unit2Unit(const Component_t& pos)
  {
    SUPERCELL_BOUNDARY_LIMITS(T);
    register T x(pos[0]),y(pos[1]),z(pos[2]);
    THREE_DIM_BOUNDARY_BLOCK(x,y,z,epsilon,plus_one);
    return TinyVector<T,3>(x,y,z);
  }

  static inline Component_t Cart2Unit(const Component_t& pos, const Transformer_t& G)
  {
    SUPERCELL_BOUNDARY_LIMITS(T);
    register T x(pos[0]*G[0]+pos[1]*G[3]+pos[2]*G[6]),
             y(pos[0]*G[1]+pos[1]*G[4]+pos[2]*G[7]),
             z(pos[0]*G[2]+pos[1]*G[5]+pos[2]*G[8]);
    THREE_DIM_BOUNDARY_BLOCK(x,y,z,epsilon,plus_one);
    return Component_t(x,y,z);
  }

  static inline Component_t Unit2Cart(const Component_t& pos, const Transformer_t& R)
  {
    SUPERCELL_BOUNDARY_LIMITS(T);
    register T x(pos[0]),y(pos[1]),z(pos[2]);
    THREE_DIM_BOUNDARY_BLOCK(x,y,z,epsilon,plus_one);
    return Component_t(x*R[0]+y*R[3]+z*R[6], x*R[1]+y*R[4]+z*R[7], x*R[2]+y*R[5]+z*R[8]);
  }

  static inline Component_t Cart2Cart(const Component_t& pos, const Transformer_t& G, const Transformer_t& R)
  {
    SUPERCELL_BOUNDARY_LIMITS(T);
    register T x(pos[0]*G[0]+pos[1]*G[3]+pos[2]*G[6]),
             y(pos[0]*G[1]+pos[1]*G[4]+pos[2]*G[7]),
             z(pos[0]*G[2]+pos[1]*G[5]+pos[2]*G[8]);
    THREE_DIM_BOUNDARY_BLOCK(x,y,z,epsilon,plus_one);
    return Component_t(x*R[0]+y*R[3]+z*R[6], x*R[1]+y*R[4]+z*R[7], x*R[2]+y*R[5]+z*R[8]);
  }
};

///////////////////////////////////////////////////
////specialization for 2D
///////////////////////////////////////////////////
template<class T>
struct ApplyBConds<ParticleAttrib<TinyVector<T,2> >, Tensor<T,2>, 2,true>
{

  typedef ParticleAttrib<TinyVector<T,2> > Array_t;
  typedef typename Array_t::Type_t         Component_t;
  typedef Tensor<T,2>                      Transformer_t;

  /** Apply boundary condition on the LatticeUnit vectors and return LatticeUnit vectors
   * @param pin input array in the LatticeUnit Unit
   * @param pout output array
   * @param first starting index
   * @param last ending index
   *
   * Move the components to [0,1)
   */
  inline static void
  Unit2Unit(const Array_t& pin, Array_t& pout, int first, int last)
  {
    SUPERCELL_BOUNDARY_LIMITS(T);
#pragma ivdep
    for(int i=first; i<last; i++)
    {
      T x(pin[i][0]),y(pin[i][1]);
      TWO_DIM_BOUNDARY_BLOCK(x,y,epsilon,plus_one);
      pout[i][0]=x;
      pout[i][1]=y;
    }
  }

  /** Apply boundary condition on the LatticeUnit vectors and return Cartesian vectors
   * @param pin input array in the LatticeUnit
   * @param R transformation matrix from LatticeUnit to CartesianUnit
   * @param pout input array in the CartesianUnit
   * @param first starting index
   * @param last ending index
   *
   * pout = dot(applybconds(pin),R)
   * - applybconds move all to [0,1)
   * - dot(X,R) convert to CartesianUnit
   */
  inline static void
  Unit2Cart(const Array_t& pin, const Transformer_t& R, Array_t& pout, int first, int last)
  {
    SUPERCELL_BOUNDARY_LIMITS(T);
    register T r00=R[0], r11=R[3];
#pragma ivdep
    for(int i=first; i<last; i++)
    {
      T x(pin[i][0]),y(pin[i][1]);
      TWO_DIM_BOUNDARY_BLOCK(x,y,epsilon,plus_one);
      pout[i][0]=r00*x;
      pout[i][1]=r11*y;
    }
  }


  /** Apply boundary condition on the Cartesin vectors and return LatticeUnit vectors
   * @param pin input array in the Cartesian Unit
   * @param G transformation matrix from CartesianUnit to LatticeUnit
   * @param pout input array in the LatticeUnit
   * @param first starting index
   * @param last ending index
   *
   * pout = applybconds(dot(pin,G))
   * - dot(pin,G) convert to LatticeUnit
   * - applybconds move all to [0,1)
   */
  inline static void
  Cart2Unit(const Array_t& pin, const Transformer_t& G, Array_t& pout, int first, int last)
  {
    SUPERCELL_BOUNDARY_LIMITS(T);
    register T g00=G[0], g11=G[3];
#pragma ivdep
    for(int i=first; i<last; i++)
    {
      T x(pin[i][0]*g00),y(pin[i][1]*g11);
      TWO_DIM_BOUNDARY_BLOCK(x,y,epsilon,plus_one);
      pout[i][0]=x;
      pout[i][1]=y;
    }
  }

  /** Apply boundary condition on the Cartesin vectors and return Cartesian vectors
   * @param pin input array in the Cartesian Unit
   * @param G transformation matrix from CartesianUnit to LatticeUnit
   * @param R transformation matrix from LatticeUnit to CartesianUnit
   * @param pout input array in the Cartesian Unit
   * @param first starting index
   * @param last ending index
   *
   * pout = dot(applybconds(dot(pin,G)),R)
   * - dot(pin,G) convert to LatticeUnit
   * - applybconds move all to [0,1)
   * - dot(X,R) convert to CartesianUnit
   */
  inline static void
  Cart2Cart(const Array_t& pin, const Transformer_t& G, const Transformer_t& R, Array_t& pout,
            int first, int last)
  {
    SUPERCELL_BOUNDARY_LIMITS(T);
    register T g00=G[0], g11=G[3];
    register T r00=R[0], r11=R[3];
#pragma ivdep
    for(int i=first; i<last; i++)
    {
      T x(pin[i][0]*g00),y(pin[i][1]*g11);
      TWO_DIM_BOUNDARY_BLOCK(x,y,epsilon,plus_one);
      pout[i][0]=r00*x;
      pout[i][1]=r11*y;
    }
  }

  /** Apply boundary condition on the LatticeUnit vectors and return LatticeUnit vectors
   * @param pinout input/output array in the LatticeUnit
   * @param first starting index
   * @param last ending index
   *
   * Move the components to [0,1)
   */
  inline static void
  Unit2Unit(Array_t& pinout, int first, int last)
  {
    SUPERCELL_BOUNDARY_LIMITS(T);
#pragma ivdep
    for(int i=first; i<last; i++)
    {
      T x(pinout[i][0]),y(pinout[i][1]);
      TWO_DIM_BOUNDARY_BLOCK(x,y,epsilon,plus_one);
      pinout[i][0]=x;
      pinout[i][1]=y;
    }
  }

  /** Apply boundary condition on the Cartesian vectors and return Cartesian vectors
   * @param pinout input/output array in the LatticeUnit
   * @param G transformation matrix from CartesianUnit to LatticeUnit
   * @param R transformation matrix from LatticeUnit to CartesianUnit
   * @param first starting index
   * @param last ending index
   *
   * pinout <- dot(applybconds(dot(pinout,G)),R)
   * - dot(pin,G) convert to LatticeUnit
   * - applybconds move all to [0,1)
   * - dot(X,R) convert to CartesianUnit
   */
  inline static void
  Cart2Cart(Array_t& pinout, const Transformer_t& G, const Transformer_t& R, int first, int last)
  {
    SUPERCELL_BOUNDARY_LIMITS(T);
    register T g00=G[0], g11=G[3];
    register T r00=R[0], r11=R[3];
#pragma ivdep
    for(int i=first; i<last; i++)
    {
      T x(pinout[i][0]*g00),y(pinout[i][1]*g11);
      TWO_DIM_BOUNDARY_BLOCK(x,y,epsilon,plus_one);
      pinout[i][0]=r00*x;
      pinout[i][1]=r11*y;
    }
  }

  static inline Component_t Unit2Unit(const Component_t& pos)
  {
    SUPERCELL_BOUNDARY_LIMITS(T);
    register T x(pos[0]),y(pos[1]);
    TWO_DIM_BOUNDARY_BLOCK(x,y,epsilon,plus_one);
    return Component_t(x,y);
  }

  static inline Component_t Cart2Unit(const Component_t& pos, const Transformer_t& G)
  {
    SUPERCELL_BOUNDARY_LIMITS(T);
    register T x(pos[0]*G[0]), y(pos[1]*G[3]);
    TWO_DIM_BOUNDARY_BLOCK(x,y,epsilon,plus_one);
    return Component_t(x,y);
  }

  static inline Component_t Unit2Cart(const Component_t& pos, const Transformer_t& R)
  {
    SUPERCELL_BOUNDARY_LIMITS(T);
    register T x(pos[0]),y(pos[1]);
    TWO_DIM_BOUNDARY_BLOCK(x,y,epsilon,plus_one);
    return Component_t(x*R[0],y*R[3]);
  }

  static inline Component_t Cart2Cart(const Component_t& pos, const Transformer_t& G, const Transformer_t& R)
  {
    SUPERCELL_BOUNDARY_LIMITS(T);
    register T x(pos[0]*G[0]),y(pos[1]*G[3]);
    TWO_DIM_BOUNDARY_BLOCK(x,y,epsilon,plus_one);
    return Component_t(x*R[0],y*R[3]);
  }
};

template<class T>
struct ApplyBConds<ParticleAttrib<TinyVector<T,2> >, Tensor<T,2>, 2,false>
{

  typedef ParticleAttrib<TinyVector<T,2> > Array_t;
  typedef typename Array_t::Type_t         Component_t;
  typedef Tensor<T,2>                      Transformer_t;

  inline static void
  Cart2Cart(const Array_t& pin, const Transformer_t& G, const Transformer_t& R, Array_t& pout, int first, int last)
  {
    SUPERCELL_BOUNDARY_LIMITS(T);
    register T g00=G[0], g01=G[1],g10=G[2],g11=G[3];
    register T r00=R[0], r01=R[1],r10=R[2],r11=R[3];
#pragma ivdep
    for(int i=first; i<last; i++)
    {
      T x(pin[i][0]*g00+pin[i][1]*g10);
      T y(pin[i][0]*g01+pin[i][1]*g11);
      TWO_DIM_BOUNDARY_BLOCK(x,y,epsilon,plus_one);
      pout[i][0]=x*r00+y*r10;
      pout[i][1]=x*r01+y*r11;
    }
  }

  inline static void
  Cart2Unit(const Array_t& pin, const Transformer_t& G, Array_t& pout, int first, int last)
  {
    SUPERCELL_BOUNDARY_LIMITS(T);
    register T g00=G[0], g01=G[1],g10=G[2],g11=G[3];
#pragma ivdep
    for(int i=first; i<last; i++)
    {
      T x(pin[i][0]*g00+pin[i][1]*g10);
      T y(pin[i][0]*g01+pin[i][1]*g11);
      TWO_DIM_BOUNDARY_BLOCK(x,y,epsilon,plus_one);
      pout[i][0]=x;
      pout[i][1]=y;
    }
  }

  inline static void
  Unit2Cart(const Array_t& pin, const Transformer_t& R, Array_t& pout, int first, int last)
  {
    SUPERCELL_BOUNDARY_LIMITS(T);
    register T r00=R[0], r01=R[1],r10=R[2],r11=R[3];
#pragma ivdep
    for(int i=first; i<last; i++)
    {
      T x(pin[i][0]), y(pin[i][1]);
      TWO_DIM_BOUNDARY_BLOCK(x,y,epsilon,plus_one);
      pout[i][0]=x*r00+y*r10;
      pout[i][1]=x*r01+y*r11;
    }
  }

  inline static void
  Unit2Unit(const Array_t& pin, Array_t& pout, int first, int last)
  {
    SUPERCELL_BOUNDARY_LIMITS(T);
#pragma ivdep
    for(int i=first; i<last; i++)
    {
      T x(pin[i][0]),y(pin[i][1]);
      TWO_DIM_BOUNDARY_BLOCK(x,y,epsilon,plus_one);
      pout[i][0]=x;
      pout[i][1]=y;
    }
  }

  inline static void
  Unit2Unit(Array_t& pinout, int first, int last)
  {
    SUPERCELL_BOUNDARY_LIMITS(T);
#pragma ivdep
    for(int i=first; i<last; i++)
    {
      T x(pinout[i][0]),y(pinout[i][1]);
      TWO_DIM_BOUNDARY_BLOCK(x,y,epsilon,plus_one);
      pinout[i][0]=x;
      pinout[i][1]=y;
    }
  }

  inline static void
  Cart2Cart(Array_t& pinout, const Transformer_t& G, const Transformer_t& R, int first, int last)
  {
    SUPERCELL_BOUNDARY_LIMITS(T);
    register T g00=G[0], g01=G[1],g10=G[2],g11=G[3];
    register T r00=R[0], r01=R[1],r10=R[2],r11=R[3];
#pragma ivdep
    for(int i=first; i<last; i++)
    {
      T x(pinout[i][0]*g00+pinout[i][1]*g10);
      T y(pinout[i][0]*g01+pinout[i][1]*g11);
      TWO_DIM_BOUNDARY_BLOCK(x,y,epsilon,plus_one);
      pinout[i][0]=x*r00+y*r10;
      pinout[i][1]=x*r01+y*r11;
    }
  }

  static inline Component_t Unit2Unit(const Component_t& pos)
  {
    SUPERCELL_BOUNDARY_LIMITS(T);
    register T x(pos[0]),y(pos[1]);
    TWO_DIM_BOUNDARY_BLOCK(x,y,epsilon,plus_one);
    return Component_t(x,y);
  }

  static inline Component_t Cart2Unit(const Component_t& pos, const Transformer_t& G)
  {
    SUPERCELL_BOUNDARY_LIMITS(T);
    register T x(pos[0]*G[0]+pos[1]*G[2]),
             y(pos[0]*G[1]+pos[1]*G[3]);
    TWO_DIM_BOUNDARY_BLOCK(x,y,epsilon,plus_one);
    return Component_t(x,y);
  }

  static inline Component_t Unit2Cart(const Component_t& pos, const Transformer_t& R)
  {
    SUPERCELL_BOUNDARY_LIMITS(T);
    register T x(pos[0]),y(pos[1]);
    TWO_DIM_BOUNDARY_BLOCK(x,y,epsilon,plus_one);
    return Component_t(x*R[0]+y*R[2], x*R[1]+y*R[3]);
  }

  static inline Component_t Cart2Cart(const Component_t& pos, const Transformer_t& G, const Transformer_t& R)
  {
    SUPERCELL_BOUNDARY_LIMITS(T);
    register T x(pos[0]*G[0]+pos[1]*G[2]),
             y(pos[0]*G[1]+pos[1]*G[3]);
    TWO_DIM_BOUNDARY_BLOCK(x,y,epsilon,plus_one);
    return Component_t(x*R[0]+y*R[2], x*R[1]+y*R[3]);
  }
};

/** inout[i]=inout[i]-floor(inout[i])
 *
 * See simd/vmath.h and should be specialized for vector libraries, e.g., INTEL vml, IBM massv
 */
template<typename T, unsigned D>
inline void put2box(ParticleAttrib<TinyVector<T,D> >& inout)
{
  simd::remainder(&(inout[0][0]),inout.size()*D);
}

/** out[i]=in[i]-floor(in[i])
 */
template<typename T, unsigned D>
inline void put2box(const ParticleAttrib<TinyVector<T,D> >& in
                    , ParticleAttrib<TinyVector<T,D> >& out)
{
  simd::remainder(&(in[0][0]),&(out[0][0]),in.size()*D);
}
}
#endif

