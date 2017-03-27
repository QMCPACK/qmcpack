//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////
    
    
/** @file Bspline3DSet.cpp
 * @brief Implement derived classes from Bspline3DBase
 */
#include "QMCWaveFunctions/Bspline3DSet.h"

namespace qmcplusplus
{

////////////////////////////////////////////////////////////
//Implementation of Bspline3DSet_Ortho
////////////////////////////////////////////////////////////
SPOSetBase* Bspline3DSet_Ortho::makeClone() const
{
  return  new Bspline3DSet_Ortho(*this);
}

void Bspline3DSet_Ortho::evaluate(const ParticleSet& e, int iat, ValueVector_t& vals)
{
  if(bKnots.Find(e.R[iat][0],e.R[iat][1],e.R[iat][2]))
  {
#pragma ivdep
    for(int j=0; j<NumOrbitals; j++)
      vals[j]=bKnots.evaluate(*P[j]);
  }
  else
  {
    vals=0.0;
  }
}

void
Bspline3DSet_Ortho::evaluate(const ParticleSet& e, int iat,
                             ValueVector_t& vals, GradVector_t& grads, ValueVector_t& laps)
{
  if(bKnots.FindAll(e.R[iat][0],e.R[iat][1],e.R[iat][2]))
  {
#pragma ivdep
    for(int j=0; j<NumOrbitals; j++)
      vals[j]=bKnots.evaluate(*P[j],grads[j],laps[j]);
  }
  else
  {
    vals=0.0;
    grads=0.0;
    laps=0.0;
  }
}

void
Bspline3DSet_Ortho::evaluate_notranspose(const ParticleSet& e, int first, int last,
    ValueMatrix_t& vals, GradMatrix_t& grads, ValueMatrix_t& laps)
{
  for(int iat=first,i=0; iat<last; iat++,i++)
  {
    if(bKnots.FindAll(e.R[iat][0],e.R[iat][1],e.R[iat][2]))
    {
#pragma ivdep
      for(int j=0; j<OrbitalSetSize; j++)
        vals(i,j)=bKnots.evaluate(*P[j],grads(i,j),laps(i,j));
    }
    else
    {
      for(int j=0; j<OrbitalSetSize; j++)
      {
        vals(i,j)=0.0;
        grads(i,j)=0.0;
        laps(i,j)=0.0;
      }
    }
  }
}


////////////////////////////////////////////////////////////
//Implementation of Bspline3DSet_Gen
////////////////////////////////////////////////////////////
SPOSetBase* Bspline3DSet_Gen::makeClone() const
{
  return new Bspline3DSet_Gen(*this);
}
void Bspline3DSet_Gen::evaluate(const ParticleSet& e, int iat, ValueVector_t& vals)
{
  PosType ru(Lattice.toUnit(e.R[iat]));
  if(bKnots.Find(ru[0],ru[1],ru[2]))
    for(int j=0; j<OrbitalSetSize; j++)
      vals[j]=bKnots.evaluate(*P[j]);
  else
    vals=0.0;
}

void
Bspline3DSet_Gen::evaluate(const ParticleSet& e, int iat,
                           ValueVector_t& vals, GradVector_t& grads, ValueVector_t& laps)
{
  PosType ru(Lattice.toUnit(e.R[iat]));
  if(bKnots.FindAll(ru[0],ru[1],ru[2]))
  {
    TinyVector<ValueType,3> gu;
    Tensor<ValueType,3> hess;
#pragma ivdep
    for(int j=0; j<OrbitalSetSize; j++)
    {
      vals[j]=bKnots.evaluate(*P[j],gu,hess);
      grads[j]=dot(Lattice.G,gu);
      laps[j]=trace(hess,GGt);
    }
  }
  else
  {
    vals=0.0;
    grads=0.0;
    laps=0.0;
  }
}

void
Bspline3DSet_Gen::evaluate_notranspose(const ParticleSet& e, int first, int last,
                                       ValueMatrix_t& vals, GradMatrix_t& grads, ValueMatrix_t& laps)
{
  for(int iat=first,i=0; iat<last; iat++,i++)
  {
    PosType ru(Lattice.toUnit(e.R[iat]));
    if(bKnots.FindAll(ru[0],ru[1],ru[2]))
    {
      TinyVector<ValueType,3> gu;
      Tensor<ValueType,3> hess;
#pragma ivdep
      for(int j=0; j<OrbitalSetSize; j++)
      {
        vals(i,j)=bKnots.evaluate(*P[j],gu,hess);
        grads(i,j)=dot(Lattice.G,gu);
        laps(i,j)=trace(hess,GGt);
      }
    }
    else
    {
      for(int j=0; j<OrbitalSetSize; j++)
      {
        vals(i,j)=0.0;
        grads(i,j)=0.0;
        laps(i,j)=0.0;
      }
      //for(int j=0; j<OrbitalSetSize; j++) vals(j,i)=0.0;
      //std::copy(grads[i],grads[i]+OrbitalSetSize,0.0);
      //std::copy(laps[i],laps[i]+OrbitalSetSize,0.0);
    }
  }
}


////////////////////////////////////////////////////////////
//Implementation of Bspline3DSet_Ortho_Trunc
////////////////////////////////////////////////////////////
SPOSetBase* Bspline3DSet_Ortho_Trunc::makeClone() const
{
  return new Bspline3DSet_Ortho_Trunc(*this);
}
void Bspline3DSet_Ortho_Trunc::evaluate(const ParticleSet& e, int iat, ValueVector_t& vals)
{
  PosType r(e.R[iat]);
  bKnots.Find(r[0],r[1],r[2]);
#pragma ivdep
  for(int j=0; j<Centers.size(); j++)
  {
    if(bKnots.getSep2(r[0]-Centers[j][0],r[1]-Centers[j][1],r[2]-Centers[j][2])>Rcut2)
      vals[j]=0.0;//numeric_limits<T>::epsilon();
    else
      vals[j]=bKnots.evaluate(*P[j]);
  }
}

void Bspline3DSet_Ortho_Trunc::evaluate(const ParticleSet& e, int iat,
                                        ValueVector_t& vals, GradVector_t& grads, ValueVector_t& laps)
{
  PosType r(e.R[iat]);
  bKnots.FindAll(r[0],r[1],r[2]);
#pragma ivdep
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

void Bspline3DSet_Ortho_Trunc::evaluate_notranspose(const ParticleSet& e, int first, int last,
    ValueMatrix_t& vals, GradMatrix_t& grads, ValueMatrix_t& laps)
{
  for(int iat=first,i=0; iat<last; iat++,i++)
  {
    PosType r(e.R[iat]);
    bKnots.FindAll(r[0],r[1],r[2]);
#pragma ivdep
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
        vals(i,j)=bKnots.evaluate(*P[j],grads(i,j),laps(i,j));
      }
    }
  }
}


////////////////////////////////////////////////////////////
//Implementation of Bspline3DSet_Gen_Trunc
////////////////////////////////////////////////////////////
SPOSetBase* Bspline3DSet_Gen_Trunc::makeClone() const
{
  return new Bspline3DSet_Gen_Trunc(*this);
}

void Bspline3DSet_Gen_Trunc::evaluate(const ParticleSet& e, int iat, ValueVector_t& vals)
{
  PosType r(e.R[iat]);
  PosType ru(Lattice.toUnit(r));
  bKnots.Find(ru[0],ru[1],ru[2]);
#pragma ivdep
  for(int j=0; j<Centers.size(); j++)
  {
    if(bKnots.getSep2(r[0]-Centers[j][0],r[1]-Centers[j][1],r[2]-Centers[j][2])>Rcut2)
      vals[j]=0.0;//numeric_limits<T>::epsilon();
    else
      vals[j]=bKnots.evaluate(*P[j]);
  }
}

void Bspline3DSet_Gen_Trunc::evaluate(const ParticleSet& e, int iat,
                                      ValueVector_t& vals, GradVector_t& grads, ValueVector_t& laps)
{
  PosType r(e.R[iat]);
  PosType ru(Lattice.toUnit(r));
  bKnots.FindAll(ru[0],ru[1],ru[2]);
  TinyVector<ValueType,3> gu;
  Tensor<ValueType,3> hess;
#pragma ivdep
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
      //vals[j]=bKnots.evaluate(*P[j],grads[j],laps[j]);
    }
  }
}


void Bspline3DSet_Gen_Trunc::evaluate_notranspose(const ParticleSet& e, int first, int last,
    ValueMatrix_t& vals, GradMatrix_t& grads, ValueMatrix_t& laps)
{
  for(int iat=first,i=0; iat<last; iat++,i++)
  {
    PosType r(e.R[iat]);
    PosType ru(Lattice.toUnit(r));
    bKnots.FindAll(ru[0],ru[1],ru[2]);
    TinyVector<ValueType,3> gu;
    Tensor<ValueType,3> hess;
#pragma ivdep
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

#if defined(QMC_COMPLEX)
////////////////////////////////////////////////////////////
//Implementation of Bspline3DSet_Twist
////////////////////////////////////////////////////////////
SPOSetBase* Bspline3DSet_Twist::makeClone() const
{
  return new Bspline3DSet_Twist(*this);
}

void Bspline3DSet_Twist::evaluate(const ParticleSet& e, int iat, ValueVector_t& vals)
{
  PosType r(e.R[iat]);
  PosType ru(Lattice.toUnit(r));
  if(bKnots.Find(ru[0],ru[1],ru[2]))
  {
    RealType phi(dot(TwistAngle,r));
    ValueType phase(std::cos(phi),std::sin(phi));
#pragma ivdep
    for(int j=0; j <OrbitalSetSize; j++)
      vals[j]=phase*bKnots.evaluate(*P[j]);
  }
  else
    vals=0.0;
}

void
Bspline3DSet_Twist::evaluate(const ParticleSet& e, int iat,
                             ValueVector_t& vals, GradVector_t& grads, ValueVector_t& laps)
{
  PosType r(e.R[iat]);
  PosType ru(Lattice.toUnit(r));
  if(bKnots.FindAll(ru[0],ru[1],ru[2]))
  {
    RealType phi(dot(TwistAngle,r));
    RealType c=std::cos(phi),s=std::sin(phi);
    ValueType phase(c,s);
    //ik e^{i{\bf k}\cdot {\bf r}}
    GradType dk(ValueType(-TwistAngle[0]*s,TwistAngle[0]*c),
                ValueType(-TwistAngle[1]*s,TwistAngle[1]*c),
                ValueType(-TwistAngle[2]*s,TwistAngle[2]*c));
    TinyVector<ValueType,3> gu;
    Tensor<ValueType,3> hess;
#pragma ivdep
    for(int j=0; j<OrbitalSetSize; j++)
    {
      ValueType v= bKnots.evaluate(*P[j],gu,hess);
      GradType g= dot(Lattice.G,gu);
      ValueType l=trace(hess,GGt);
      vals[j]=phase*v;
      grads[j]=v*dk+phase*g;
      laps[j]=phase*(mK2*v+l)+2.0*dot(dk,g);
    }
  }
  else
  {
    vals=0.0;
    grads=0.0;
    laps=0.0;
  }
}

void
Bspline3DSet_Twist::evaluate_notranspose(const ParticleSet& e, int first, int last,
    ValueMatrix_t& vals, GradMatrix_t& grads, ValueMatrix_t& laps)
{
  for(int iat=first,i=0; iat<last; iat++,i++)
  {
    PosType r(e.R[iat]);
    PosType ru(Lattice.toUnit(r));
    if(bKnots.FindAll(ru[0],ru[1],ru[2]))
    {
      RealType phi(dot(TwistAngle,r));
      RealType c=std::cos(phi),s=std::sin(phi);
      ValueType phase(c,s);
      GradType dk(ValueType(-TwistAngle[0]*s,TwistAngle[0]*c),
                  ValueType(-TwistAngle[1]*s,TwistAngle[1]*c),
                  ValueType(-TwistAngle[2]*s,TwistAngle[2]*c));
      TinyVector<ValueType,3> gu;
      Tensor<ValueType,3> hess;
#pragma ivdep
      for(int j=0; j<OrbitalSetSize; j++)
      {
        ValueType v=bKnots.evaluate(*P[j],gu,hess);
        GradType g=dot(Lattice.G,gu);
        ValueType l=trace(hess,GGt);
        vals(i,j)=phase*v;
        grads(i,j)=v*dk+phase*g;
        laps(i,j)=phase*(mK2*v+l)+2.0*dot(dk,g);
      }
    }
    else
    {
      for(int j=0; j<OrbitalSetSize; j++)
      {
        vals(i,j)=0.0;
        grads(i,j)=0.0;
        laps(i,j)=0.0;
      }
    }
  }
}
#endif
}
