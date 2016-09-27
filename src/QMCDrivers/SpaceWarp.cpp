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
    
    


#include "QMCDrivers/SpaceWarp.h"
#include "QMCDrivers/SinglePtclWarp.h"

namespace qmcplusplus
{

void SpaceWarp::update_one_ptcl_Jacob(int iel)
{
  for(int ipsi=0; ipsi< npsi; ipsi++)
    one_ptcl_Jacob[ipsi][iel]=WarpVector[iel]->Jacobian[ipsi];
}

//void SpaceWarp::initialize(std::vector<DistanceTableData*>& dtList){
void SpaceWarp::initialize(std::vector<ParticleSet*>& ionSets, DistanceTableData* dtprime)
{
  dtPrimary=dtprime;
  npsi=ionSets.size();
  std::cout << dtPrimary->VisitorIndex << std::endl;
  nptcl= dtPrimary->size(DistanceTableData::VisitorIndex);
  ncenter=dtPrimary->centers();
  //resize and Initialize nuclear displacement
  Delta.resize(npsi);
  for(int ipsi=0; ipsi<npsi; ipsi++)
  {
    Delta[ipsi].resize(ncenter);
  }
  for (int iptcl=0; iptcl < nptcl; iptcl++)
    WarpVector.push_back(new SinglePtclWarp(ncenter,npsi));
  //Compute displacements of ions for each ionic configuration
  for(int ipsi=0; ipsi<npsi; ipsi++)
  {
    for(int iat=0; iat<ncenter; iat++)
      //(Delta[ipsi])[iat]= dtList[ipsi]->origin().R[iat] - dtPrimary->origin().R[iat];
      (Delta[ipsi])[iat]= ionSets[ipsi]->R[iat] - dtPrimary->origin().R[iat];
  }
  //auxiliary vectors
  r.resize(ncenter);
  rinv.resize(ncenter);
  dr.resize(ncenter);
}

void SpaceWarp::warp_one(int iel,bool require_register)
{
  if(require_register)
  {
    for(int iat=0; iat<ncenter; iat++)
    {
      r[iat]=dtPrimary->Temp[iat].r1;
      rinv[iat]=dtPrimary->Temp[iat].rinv1;
      dr[iat]=dtPrimary->Temp[iat].dr1;
    }
  }
  else
  {
    for(int iat=0; iat<ncenter; iat++)
    {
      int nn=iat*nptcl+iel;
      r[iat]=dtPrimary->r(nn);
      rinv[iat]=dtPrimary->rinv(nn);
      dr[iat]=dtPrimary->dr(nn);
    }
  }
  WarpVector[iel]->update_warp_disp(r,rinv,dr);
}

QMCTraits::PosType SpaceWarp::get_displacement(int iel,int ipsi)
{
  return WarpVector[iel]->get_displacement(ipsi,Delta[ipsi]);
}

QMCTraits::RealType SpaceWarp::get_Jacobian(int iel, int ipsi)
{
  WarpVector[iel]->update_warp_jacob();
  return WarpVector[iel]->get_Jacobian(ipsi,Delta[ipsi]);
}

QMCTraits::TensorType SpaceWarp::get_Jacob_matrix(int iel, int ipsi)
{
  return WarpVector[iel]->Jacob_matrix[ipsi];
}

QMCTraits::TensorType SpaceWarp::get_Jacob_cofactor(int iel, int ipsi)
{
  return WarpVector[iel]->Jacob_cofactor[ipsi];
}

QMCTraits::PosType SpaceWarp::get_grad_ln_Jacob(int iel, int ipsi)
{
  WarpVector[iel]->update_warp_grad();
  return WarpVector[iel]->get_grad_ln_Jacob(ipsi,Delta[ipsi]);
}

QMCTraits::PosType SpaceWarp::get_grad_ln_Jacob_num(int iel, int ipsi)
{
  PosType grad_ln_Jacob_num;
  RealType epsilon=1.0e-4;
  RealType direction[2],logJacob[2];
  direction[0]=-1.e0;
  direction[1]=1.e0;
  for(int igamma=0; igamma<3; igamma++)
  {
    for(int id=0; id<2; id++)
    {
      for(int iat=0; iat<ncenter; iat++)
      {
        dr[iat]=0.e0;
        (dr[iat])[igamma]=direction[id]*epsilon;
        int nn=iat*nptcl+iel;
        dr[iat]+=dtPrimary->dr(nn);
        r[iat]=0.e0;
        for(int ialpha=0; ialpha<3; ialpha++)
          r[iat]+=pow((dr[iat])[ialpha],2);
        r[iat]=sqrt(r[iat]);
        rinv[iat]=1.e0/r[iat];
      }
      WarpVector[iel]->update_warp_disp(r,rinv,dr);
      WarpVector[iel]->get_displacement(ipsi,Delta[ipsi]);
      WarpVector[iel]->update_warp_jacob();
      logJacob[id]=log(WarpVector[iel]->get_Jacobian(ipsi,Delta[ipsi]));
    }
    grad_ln_Jacob_num[igamma]=0.5*(logJacob[1]-logJacob[0])/epsilon;
  }
  return grad_ln_Jacob_num;
}

void SpaceWarp::registerData(std::vector<ParticleSet*>& plist, PooledData<RealType>& buf)
{
  if(PtclRefs.empty())
  {
    SizeOfR=plist[0]->getTotalNum()*DIM;
    for(int ipsi=0; ipsi<npsi; ipsi++)
    {
      PtclRefs.push_back(plist[ipsi]);
    }
  }
  for(int ipsi=0; ipsi<npsi; ipsi++)
  {
    RealType* first=&(plist[ipsi]->R[0][0]);
    buf.add(first,first+SizeOfR);
  }
  buf.add(one_ptcl_Jacob.begin(),one_ptcl_Jacob.end());
}

void SpaceWarp::updateBuffer(PooledData<RealType>& buf)
{
//recompute Jacobian from scratch
  copyToBuffer(buf);
}

void SpaceWarp::copyToBuffer(PooledData<RealType>& buf)
{
  std::vector<ParticleSet*>::iterator pit(PtclRefs.begin()), pit_end(PtclRefs.end());
  while(pit != pit_end)
  {
    RealType* first=&((**pit).R[0][0]);
    buf.put(first,first+SizeOfR);
    ++pit;
  }
  buf.put(one_ptcl_Jacob.begin(),one_ptcl_Jacob.end());
}

void SpaceWarp::copyFromBuffer(PooledData<RealType>& buf)
{
  std::vector<ParticleSet*>::iterator pit(PtclRefs.begin()), pit_end(PtclRefs.end());
  while(pit != pit_end)
  {
    RealType* first=&((**pit).R[0][0]);
    buf.get(first,first+SizeOfR);
    ++pit;
  }
  buf.get(one_ptcl_Jacob.begin(),one_ptcl_Jacob.end());
}
}

