//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////
    
    


#include "ParticleIO/ParticleIOUtility.h"
#include "Utilities/ProgressReportEngine.h"
namespace qmcplusplus
{

#if OHMMS_DIM ==3
void expandSuperCell(ParticleSet& ref_, Tensor<int,3>& tmat)
{
  typedef ParticleSet::SingleParticlePos_t SingleParticlePos_t;
  typedef ParticleSet::Tensor_t Tensor_t;
  Tensor<int,3> I(1,0,0,0,1,0,0,0,1);
  bool identity=true;
  int ij=0;
  while(identity&& ij<9)
  {
    identity=(I[ij]==tmat[ij]);
    ++ij;
  }
  if(identity)
    return;
  ReportEngine PRE("expandSuperCell"," ");
  app_log() << "  TileMatrix != Identity. Expanding a simulation cell for "
            << ref_.getName() << std::endl;
  {
    char buff[500];
    snprintf (buff, 500
              , "   tilematrix= %4d %4d %4d %4d %4d %4d %4d %4d %4d\n"
              , tmat[0], tmat[1], tmat[2] ,tmat[3], tmat[4], tmat[5], tmat[6], tmat[7], tmat[8]);
    app_log() << buff<< std::endl;
  }
  //convert2unit
  ref_.convert2Unit(ref_.R);
  ParticleSet::ParticleLayout_t PrimCell(ref_.Lattice);
  ref_.Lattice.set(dot(tmat,PrimCell.R));
  int natoms=ref_.getTotalNum();
  int numCopies = std::abs(det(tmat));
  ParticleSet::ParticlePos_t primPos(ref_.R);
  ParticleSet::ParticleIndex_t primTypes(ref_.GroupID);
  ref_.resize(natoms*numCopies);
  int maxCopies = 10;
  int index=0;
  //set the unit to the Cartesian
  ref_.R.InUnit=PosUnit::CartesianUnit;
  app_log() << "  Reduced coord    Cartesion coord    species.\n";
  for(int ns=0; ns<ref_.getSpeciesSet().getTotalNum(); ++ns)
  {
    for (int i0=-maxCopies; i0<=maxCopies; i0++)
      for (int i1=-maxCopies; i1<=maxCopies; i1++)
        for (int i2=-maxCopies; i2<=maxCopies; i2++)
          for (int iat=0; iat < primPos.size(); iat++)
          {
            if(primTypes[iat]!=ns)
              continue;
            //SingleParticlePos_t r     = primPos[iat];
            SingleParticlePos_t uPrim = primPos[iat];
            for (int i=0; i<3; i++)
              uPrim[i] -= std::floor(uPrim[i]);
            SingleParticlePos_t r = PrimCell.toCart(uPrim) + (double)i0*PrimCell.a(0)
                                    + (double)i1*PrimCell.a(1) + (double)i2*PrimCell.a(2);
            SingleParticlePos_t uSuper = ref_.Lattice.toUnit(r);
            if ((uSuper[0] >= -1.0e-6) && (uSuper[0] < 0.9999) &&
                (uSuper[1] >= -1.0e-6) && (uSuper[1] < 0.9999) &&
                (uSuper[2] >= -1.0e-6) && (uSuper[2] < 0.9999))
            {
              char buff[500];
              snprintf (buff, 500, "  %10.4f  %10.4f %10.4f   %12.6f %12.6f %12.6f %d\n",
                        uSuper[0], uSuper[1], uSuper[2], r[0], r[1], r[2], ns);
              app_log() << buff;
              ref_.R[index]= r;
              ref_.GroupID[index]= ns;//primTypes[iat];
              ref_.ID[index]=index;
              ref_.PCID[index]=iat;
              index++;
            }
          }
  }
  app_log() << "  Simulationcell after tiling" << std::endl;
  ref_.Lattice.print(app_log());
  app_log() << std::endl;
}

#elif OHMMS_DIM == 2
void expandSuperCell(ParticleSet& ref_, Tensor<int,2>& tmat)
{
  typedef ParticleSet::SingleParticlePos_t SingleParticlePos_t;
  typedef ParticleSet::Tensor_t Tensor_t;
  Tensor<int,2> I(1,0,0,1);
  bool identity=true;
  int ij=0;
  while(identity&& ij<4)
  {
    identity=(I[ij]==tmat[ij]);
    ++ij;
  }
  if(identity)
    return;
  //convert2unit
  ref_.convert2Unit(ref_.R);
  ParticleSet::ParticleLayout_t PrimCell(ref_.Lattice);
  ref_.Lattice.set(dot(tmat,PrimCell.R));
  int natoms=ref_.getTotalNum();
  int numCopies = std::abs(det(tmat));
  ParticleSet::ParticlePos_t primPos(ref_.R);
  ParticleSet::ParticleIndex_t primTypes(ref_.GroupID);
  ref_.resize(natoms*numCopies);
  int maxCopies = 10;
  int index=0;
  //set the unit to the Cartesian
  ref_.R.InUnit=PosUnit::CartesianUnit;
  app_log() << "  Reduced coord    Cartesion coord    species.\n";
  for(int ns=0; ns<ref_.getSpeciesSet().getTotalNum(); ++ns)
  {
    for (int i0=-maxCopies; i0<=maxCopies; i0++)
      for (int i1=-maxCopies; i1<=maxCopies; i1++)
        for (int iat=0; iat < primPos.size(); iat++)
        {
          if(primTypes[iat]!=ns)
            continue;
          //SingleParticlePos_t r     = primPos[iat];
          SingleParticlePos_t uPrim = primPos[iat];
          for (int i=0; i<2; i++)
            uPrim[i] -= std::floor(uPrim[i]);
          SingleParticlePos_t r = PrimCell.toCart(uPrim) + (double)i0*PrimCell.a(0)
                                  + (double)i1*PrimCell.a(1);
          SingleParticlePos_t uSuper = ref_.Lattice.toUnit(r);
          if ((uSuper[0] >= -1.0e-6) && (uSuper[0] < 0.9999) &&
              (uSuper[1] >= -1.0e-6) && (uSuper[1] < 0.9999) )
          {
            char buff[500];
            snprintf (buff, 500, "  %10.4f  %10.4f   %12.6f %12.6f %d\n", uSuper[0], uSuper[1], r[0], r[1],  ns);
            app_log() << buff;
            ref_.R[index]= r;
            ref_.GroupID[index]= ns;//primTypes[iat];
            index++;
          }
        }
  }
}
#else
#error "Only 2D and 3D are implemented in ParticleIO/ParticleIOUtilitcy.cpp"
#endif
}

