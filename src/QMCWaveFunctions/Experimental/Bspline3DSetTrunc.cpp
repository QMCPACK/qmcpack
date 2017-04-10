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
    
    
/** @file Bspline3DSet.cpp
 * @brief Implement derived classes from Bspline3DBase
 */
#include "QMCWaveFunctions/Bspline3DSetTrunc.h"

namespace qmcplusplus
{

////////////////////////////////////////////////////////////
//Implementation of Bspline3DSet_MLW
////////////////////////////////////////////////////////////
Bspline3DSet_MLW::Bspline3DSet_MLW()
{
  //testing only
  Lx=1.3665365600e+01;
  Ly=1.3665365600e+01;
  Lz=1.3665365600e+01;
  LxInv=1.0/Lx;
  LyInv=1.0/Ly;
  LzInv=1.0/Lz;
  LxSq=Lx*Lx;
  LySq=Ly*Ly;
  LzSq=Lz*Lz;
}

Bspline3DSet_MLW::~Bspline3DSet_MLW()
{
}

void Bspline3DSet_MLW::evaluate(const ParticleSet& e, int iat, ValueVector_t& vals)
{
  PosType r(e.R[iat]);
  for(int j=0; j<OrbitalSetSize; j++)
  {
    if(getSep2(r[0]-Centers[j][0],r[1]-Centers[j][1],r[2]-Centers[j][2])>Rcut2)
      vals[j]=0.0;
    else
    {
      PosType rtr=translate(r,j);
      bKnots.Find(rtr[0],rtr[1],rtr[2]);
      vals[j]=bKnots.evaluate(*P[j]);
    }
  }
}

void Bspline3DSet_MLW::evaluate(const ParticleSet& e, int iat,
                                ValueVector_t& vals, GradVector_t& grads, ValueVector_t& laps)
{
  PosType r(e.R[iat]);
  for(int j=0; j<OrbitalSetSize; j++)
  {
    if(getSep2(r[0]-Centers[j][0],r[1]-Centers[j][1],r[2]-Centers[j][2])>Rcut2)
    {
      vals[j]=0.0;
      grads[j]=0.0;
      laps[j]=0.0;
    }
    else
    {
      PosType rtr=translate(r,j);
      bKnots.FindAll(rtr[0],rtr[1],rtr[2]);
      vals[j]=bKnots.evaluate(*P[j],grads[j],laps[j]);
    }
  }
}

void Bspline3DSet_MLW::evaluate(const ParticleSet& e, int first, int last,
                                ValueMatrix_t& vals, GradMatrix_t& grads, ValueMatrix_t& laps)
{
  for(int iat=first,i=0; iat<last; iat++,i++)
  {
    PosType r(e.R[iat]);
    for(int j=0; j<OrbitalSetSize; j++)
    {
      //PosType dr=r-Centers[j];
      //RealType sep2=getSep2(dr[0],dr[1],dr[2]);
      //if(sep2>Rcut2)
      if(getSep2(r[0]-Centers[j][0],r[1]-Centers[j][1],r[2]-Centers[j][2])>Rcut2)
      {
        //cout << "Skip " << j << " | " << sep2 << " | " << Centers[j] << std::endl;
        vals(j,i)=0.0;
        grads(i,j)=0.0;
        laps(i,j)=0.0;
      }
      else
      {
        PosType rtr(translate(r,j));
        bKnots.FindAll(rtr[0],rtr[1],rtr[2]);
        vals(j,i)=bKnots.evaluate(*P[j],grads(i,j),laps(i,j));
        //cout << "Skip " << j << " | " << sep2 << " | " << Centers[j] << std::endl;
      }
    }
  }
}

void Bspline3DSet_MLW::evaluate_notranspose(const ParticleSet& e, int first, int last,
    ValueMatrix_t& vals, GradMatrix_t& grads, ValueMatrix_t& laps)
{
  for(int iat=first,i=0; iat<last; iat++,i++)
  {
    PosType r(e.R[iat]);
    for(int j=0; j<OrbitalSetSize; j++)
    {
      //PosType dr=r-Centers[j];
      //RealType sep2=getSep2(dr[0],dr[1],dr[2]);
      //if(sep2>Rcut2)
      if(getSep2(r[0]-Centers[j][0],r[1]-Centers[j][1],r[2]-Centers[j][2])>Rcut2)
      {
        //cout << "Skip " << j << " | " << sep2 << " | " << Centers[j] << std::endl;
        vals(i,j)=0.0;
        grads(i,j)=0.0;
        laps(i,j)=0.0;
      }
      else
      {
        PosType rtr(translate(r,j));
        bKnots.FindAll(rtr[0],rtr[1],rtr[2]);
        vals(i,j)=bKnots.evaluate(*P[j],grads(i,j),laps(i,j));
        //cout << "Skip " << j << " | " << sep2 << " | " << Centers[j] << std::endl;
      }
    }
  }
}
}
