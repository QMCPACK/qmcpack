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
/** @file Bspline3DSet.cpp
 * @brief Implement Bspline3DBase and its derived classes
 */
#include "QMCWaveFunctions/Bspline3DSet.h"

namespace qmcplusplus {

  ////////////////////////////////////////////////////////////
  //Implementation of Bspline3DSet_Ortho
  ////////////////////////////////////////////////////////////
  void Bspline3DSet_Ortho::evaluate(const PosType& r, ValueVector_t& vals) 
  {
    bKnots.Find(r[0],r[1],r[2]);
    for(int j=0; j<NumOrbitals; j++) 
      vals[j]=bKnots.evaluate(*P[j]);
  }

  void
    Bspline3DSet_Ortho::evaluate(const PosType& r, ValueVector_t& vals, 
        GradVector_t& grads, ValueVector_t& laps)
    {
      bKnots.FindAll(r[0],r[1],r[2]);
      for(int j=0; j<NumOrbitals; j++) 
        vals[j]=bKnots.evaluate(*P[j],grads[j],laps[j]);
    }

  void
    Bspline3DSet_Ortho::evaluate(const PosType& r, int i, ValueMatrix_t& vals, 
        GradMatrix_t& grads, ValueMatrix_t& laps)
    {
      bKnots.FindAll(r[0],r[1],r[2]);
      for(int j=0; j<NumOrbitals; j++) 
        vals(j,i)=bKnots.evaluate(*P[j],grads(i,j),laps(i,j));
    }

  ////////////////////////////////////////////////////////////
  //Implementation of Bspline3DSet_Gen
  ////////////////////////////////////////////////////////////
  void Bspline3DSet_Gen::evaluate(const PosType& r, ValueVector_t& vals) 
  {
    PosType ru(Lattice.toUnit(r));
    bKnots.Find(ru[0],ru[1],ru[2]);
    for(int j=0; j<NumOrbitals; j++) vals[j]=bKnots.evaluate(*P[j]);
  }

  void
    Bspline3DSet_Gen::evaluate(const PosType& r, ValueVector_t& vals, 
        GradVector_t& grads, ValueVector_t& laps)
    {
      PosType ru(Lattice.toUnit(r));
      TinyVector<value_type,3> gu;
      Tensor<value_type,3> hess;
      bKnots.FindAll(ru[0],ru[1],ru[2]);
      for(int j=0; j<NumOrbitals; j++)
      {
        vals[j]=bKnots.evaluate(*P[j],gu,hess);
        grads[j]=dot(Lattice.G,gu);
        laps[j]=trace(hess,GGt);
      }
    }

  void
    Bspline3DSet_Gen::evaluate(const PosType& r, int i, ValueMatrix_t& vals, 
        GradMatrix_t& grads, ValueMatrix_t& laps)
    {
      PosType ru(Lattice.toUnit(r));
      TinyVector<value_type,3> gu;
      Tensor<value_type,3> hess;
      bKnots.FindAll(ru[0],ru[1],ru[2]);
      for(int j=0; j<NumOrbitals; j++)
      {
        vals(j,i)=bKnots.evaluate(*P[j],gu,hess);
        grads(i,j)=dot(Lattice.G,gu);
        laps(i,j)=trace(hess,GGt);
      }
    }

  ////////////////////////////////////////////////////////////
  //Implementation of Bspline3DSet_Ortho_Trunc
  ////////////////////////////////////////////////////////////
  void Bspline3DSet_Ortho_Trunc::evaluate(const PosType& r, ValueVector_t& vals) 
  {
    bKnots.Find(r[0],r[1],r[2]);
    for(int j=0; j<NumOrbitals; j++) {
      if(bKnots.getSep2(r[0]-Centers[j][0],r[1]-Centers[j][1],r[2]-Centers[j][2])>Rcut2) 
        vals[j]=0.0;//numeric_limits<T>::epsilon();
      else
        vals[j]=bKnots.evaluate(*P[j]);
    }
  }

  void
    Bspline3DSet_Ortho_Trunc::evaluate(const PosType& r, ValueVector_t& vals, 
        GradVector_t& grads, ValueVector_t& laps)
    {
      bKnots.FindAll(r[0],r[1],r[2]);
      for(int j=0; j<NumOrbitals; j++) {
        if(bKnots.getSep2(r[0]-Centers[j][0],r[1]-Centers[j][1],r[2]-Centers[j][2])>Rcut2) 
        {
          vals[j]=0.0;//numeric_limits<T>::epsilon();
          grads[j]=0.0;laps[j]=0.0;
        }
        else
          vals[j]=bKnots.evaluate(*P[j],grads[j],laps[j]);
      }
    }

  void
    Bspline3DSet_Ortho_Trunc::evaluate(const PosType& r, int i, ValueMatrix_t& vals, 
        GradMatrix_t& grads, ValueMatrix_t& laps)
    {
      bKnots.FindAll(r[0],r[1],r[2]);
      for(int j=0; j<NumOrbitals; j++) {
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

  ////////////////////////////////////////////////////////////
  //Implementation of Bspline3DSet_Gen_Trunc
  ////////////////////////////////////////////////////////////
  void Bspline3DSet_Gen_Trunc::evaluate(const PosType& r, ValueVector_t& vals) 
  {
    PosType ru(Lattice.toUnit(r));
    bKnots.Find(ru[0],ru[1],ru[2]);
    for(int j=0; j<NumOrbitals; j++) {
      if(bKnots.getSep2(r[0]-Centers[j][0],r[1]-Centers[j][1],r[2]-Centers[j][2])>Rcut2) 
        vals[j]=0.0;//numeric_limits<T>::epsilon();
      else
        vals[j]=bKnots.evaluate(*P[j]);
    }
  }

  void
    Bspline3DSet_Gen_Trunc::evaluate(const PosType& r, ValueVector_t& vals, 
        GradVector_t& grads, ValueVector_t& laps)
    {
      PosType ru(Lattice.toUnit(r));
      bKnots.FindAll(ru[0],ru[1],ru[2]);
      TinyVector<value_type,3> gu;
      Tensor<value_type,3> hess;
      for(int j=0; j<NumOrbitals; j++) {
        if(bKnots.getSep2(r[0]-Centers[j][0],r[1]-Centers[j][1],r[2]-Centers[j][2])>Rcut2) 
        {
          vals[j]=0.0;//numeric_limits<T>::epsilon();
          grads[j]=0.0;laps[j]=0.0;
        }
        else
        {
          vals[j]=bKnots.evaluate(*P[j],gu,hess);
          grads[j]=dot(Lattice.G,gu);
          laps[j]=trace(hess,GGt);
          //vals[j]=bKnots.evaluate(*P[j],grads[j],laps[j]);
        }
      }
    }

  void
    Bspline3DSet_Gen_Trunc::evaluate(const PosType& r, int i, ValueMatrix_t& vals, 
        GradMatrix_t& grads, ValueMatrix_t& laps)
    {
      PosType ru(Lattice.toUnit(r));
      bKnots.FindAll(ru[0],ru[1],ru[2]);
      TinyVector<value_type,3> gu;
      Tensor<value_type,3> hess;
      for(int j=0; j<NumOrbitals; j++) {
        if(bKnots.getSep2(r[0]-Centers[j][0],r[1]-Centers[j][1],r[2]-Centers[j][2])>Rcut2) 
        {
          vals(j,i)=0.0; //numeric_limits<T>::epsilon();
          grads(i,j)=0.0;
          laps(i,j)=0.0;
        }
        else
        {
          vals(j,i)=bKnots.evaluate(*P[j],gu,hess);
          grads(i,j)=dot(Lattice.G,gu);
          laps(i,j)=trace(hess,GGt);
          //vals(j,i)=bKnots.evaluate(*P[j],grads(i,j),laps(i,j));
        }
      }
    }

#if defined(QMC_COMPLEX)
  ////////////////////////////////////////////////////////////
  //Implementation of Bspline3DSet_Twist
  ////////////////////////////////////////////////////////////
  void Bspline3DSet_Twist::evaluate(const PosType& r, ValueVector_t& vals) 
  {
    PosType ru(Lattice.toUnit(r));
    bKnots.Find(ru[0],ru[1],ru[2]);
    real_type phi(dot(TwistAngle,r));
    value_type phase(std::cos(phi),std::sin(phi));
    for(int j=0; j <NumOrbitals; j++)
      vals[j]=phase*bKnots.evaluate(*P[j]);
  }

  void
    Bspline3DSet_Twist::evaluate(const PosType& r, ValueVector_t& vals, 
        GradVector_t& grads, ValueVector_t& laps)
    {
      PosType ru(Lattice.toUnit(r));
      bKnots.FindAll(ru[0],ru[1],ru[2]);

      real_type phi(dot(TwistAngle,r));
      real_type c=std::cos(phi),s=std::sin(phi);
      value_type phase(c,s);

      //ik e^{i{\bf k}\cdot {\bf r}}
      GradType dk(value_type(-TwistAngle[0]*s,TwistAngle[0]*c),
          value_type(-TwistAngle[1]*s,TwistAngle[1]*c),
          value_type(-TwistAngle[2]*s,TwistAngle[2]*c));
      TinyVector<value_type,3> gu;
      Tensor<value_type,3> hess;
      for(int j=0; j<NumOrbitals; j++) 
      {
        value_type v= bKnots.evaluate(*P[j],gu,hess);
        GradType g= dot(Lattice.G,gu);
        value_type l=trace(hess,GGt);
        vals[j]=phase*v;
        grads[j]=v*dk+phase*g;
        laps[j]=phase*(mK2*v+l)+2.0*dot(dk,g);
      }
    }

  void
    Bspline3DSet_Twist::evaluate(const PosType& r, int i, ValueMatrix_t& vals, 
        GradMatrix_t& grads, ValueMatrix_t& laps)
    {
      PosType ru(Lattice.toUnit(r));
      bKnots.FindAll(ru[0],ru[1],ru[2]);

      real_type phi(dot(TwistAngle,r));
      real_type c=std::cos(phi),s=std::sin(phi);
      value_type phase(c,s);

      GradType dk(value_type(-TwistAngle[0]*s,TwistAngle[0]*c),
          value_type(-TwistAngle[1]*s,TwistAngle[1]*c),
          value_type(-TwistAngle[2]*s,TwistAngle[2]*c));
      TinyVector<value_type,3> gu;
      Tensor<value_type,3> hess;
      for(int j=0; j<NumOrbitals; j++) 
      {
        value_type v=bKnots.evaluate(*P[j],gu,hess);
        GradType g=dot(Lattice.G,gu);
        value_type l=trace(hess,GGt);
        vals(j,i)=phase*v;
        grads(i,j)=v*dk+phase*g;
        laps(i,j)=phase*(mK2*v+l)+2.0*dot(dk,g);
      }
    }
#endif
}
/***************************************************************************
 * $RCSfile$   $Author: jnkim $
 * $Revision: 2013 $   $Date: 2007-05-22 16:47:09 -0500 (Tue, 22 May 2007) $
 * $Id: TricubicBsplineSet.h 2013 2007-05-22 21:47:09Z jnkim $
 ***************************************************************************/
