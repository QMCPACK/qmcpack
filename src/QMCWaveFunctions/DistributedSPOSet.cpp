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
    
    
#include "QMCWaveFunctions/DistributedSPOSet.h"
//this will go into cmake
#define MAX_NUM_SHARED_NODES 8
#define PSI_DIM 5

namespace qmcplusplus
{

///constructor
DistributedSPOSet::DistributedSPOSet(int norbs=0)
{
  setOrbitalSetSize(norbs);
}

DistributedSPOSet::~DistributedSPOSet()
{
  delete_iter(SendBuffer.begin(),SendBuffer.end());
  delete_iter(RecvBuffer.begin(),RecvBuffer.end());
}

void DistributedSPOSet::setCommunicator(Communicate* c)
{
  if(myComm && myComm == c)
  {
    app_log() << "  Identical communicator. Nothing to be done in ScalarEstimatorManager::setCommunicator." << std::endl;
    return;
  }
  if(c)
    myComm=c;
  else
    myComm = OHMMS::Controller;
  int ncontexts=myComm->ncontexts();
  Rnow.resize(OrbitalSetSize);
  //divide orbitals among the processors
  OrbitalCount.resize(ncontexts);
  OrbitalOffset.resize(ncontexts+1);
  FairDivideLow(OrbitalSetSize,ncontexts,OrbitalOffSet);
  for(int i=0; i<ncontexts; i++)
    OrbitalCount[i]=OrbitalOffset[i+1]-OrbitalOffset[i];
  NumRemoteNodes=ncontexts-1;
  if(RemoteNodes.empty())//cannot do it again
  {
    RemoteNodes.resize(NumRemoteNodes);
    SendBuffer.resize(ncontexts,0);
    RecvBuffer.resize(ncontexts,0);
    for(int i=0, r=0; i<=NumRemoteNodes; i++)
    {
      if(i == myNodeID)
      {
        SendBuffer[i]=new BufferType;
        RecvBuffer[i]=new BufferType;
      }
      else
      {
        RemonteNodes[r++]=i;
        SendBuffer[i]=new BufferType(OrbitalCount[i]*OrbitalSetSize*PSI_DIM);
        RecvBuffer[i]=new BufferType(OrbitalCount[i]*OrbitalSetSize*PSI_DIM);
      }
    }
  }
}

void DistributedSPOSet::setOrbitalSetSize(int norbs)
{
  if(norbs == OrbitalSetSize )
    return;
  OrbitalSetSize=norbs;
  BasisSetSize=norbs;
}

void DistributedSPOSet::resetParameters(VarRegistry<RealType>& optVariables)
{
  Phi->resetParameters(optVariables);
}

void DistributedSPOSet::resetTargetParticleSet(ParticleSet& P)
{
  Phi->resetTargetParticleSet(P);
}

void
DistributedSPOSet::evaluate(const ParticleSet& P, int iat, ValueVector_t& psi)
{
  CommunicatorTraits::mpi_request_type sendPos[MAX_NUM_SHARED_NODES];
  CommunicatorTraits::mpi_request_type recvPsi[MAX_NUM_SHARED_NODES];
  CommunicatorTraits::mpi_request_type sendPsi[MAX_NUM_SHARED_NODES];
  CommunicatorTraits::mpi_status_type statusPos[MAX_NUM_SHARED_NODES];
  CommunicatorTraits::mpi_status_type statusPsi[MAX_NUM_SHARED_NODES];
  PosType pos(P[iat]);
  //send the current position and start recv
  for(int p=0; p<NumRemoteNodes; p++)
  {
    int target=RemoteNodes[p];
    MPI_Isend(pos.begin(), OHMMS_DIM, MPI_DOUBLE, target, PosTag, myComm->getMPI(), &(sendPos[p]));
    MPI_Irecv(RecvBuffer[target]->data(), OrbitalCount[target], target, PsiTag, myComm->getMPI(),&(recvPsi[p]));
  }
  ValueVector_t psiL(OrbitalCount[myNodeID]);
  //do the local calculation
  Phi->evaluate(pos,psiL);
  copy(psiL.begin(),psiL.end(),psi.begin()+OrbitalOffset[myNodeID]);
  //can make waitsome
  int err=MPI_Waitall(NumRemoteNodes, sendPos,statusPos);
  //do the calculation with the positions recv and send back the orbitals
  for(int p=0; p<NumRemoteNodes; p++)
  {
    int target=RemoteNodes[p];
    //matching recv with MPI_Isend
    MPI_recv(Rnow[p].data(), OHMMS_DIM, MPI_DOUBLE, target, PosTag, myNodeID,myComm->getMPI(),&(statusPos[p]));
    //calculate with a position
    Phi->evaluate(Rnow[p],psiL);
    //pack the message
    SendBuffer[p]->rewind();
    SendBuffer[p]->put(psiL.begin(),psiL.end());
    //send the wavefunction matching ready-send with MPI_Irecv
    MPI_Irsend(SendBuffer[p]->data(), OrbitalCount[myNodeID], MPI_DOUBLE, target, PsiTag, myComm->getMPI(), &(sendPsi[p]));
  }
  int nm=NumRemoteNodes;
  while(nm)
  {
    int count =0;
    int err=MPI_Testsome(NumRemoteNodes,recvPsi,&count,statusPsi);
    for(int m=0; m<count; m++)
    {
      int source=statusPsi[m].MPI_SOURCE;
      RecvBuffer[source]->rewind();
      for(int t=OrbitalOffset[source]; t<OrbitalOffset[source+1]; t++)
        RecvBuffer[source]->get(psi[t]);
    }
    nm-=count;
  }
  err=MPI_Waitall(NumRemoteNodes, sendPsi,statusPsi);
}

void
DistributedSPOSet::evaluate(const ParticleSet& P, int iat,
                            ValueVector_t& psi, GradVector_t& dpsi, ValueVector_t& d2psi)
{
  CommunicatorTraits::mpi_request_type sendPos[MAX_NUM_SHARED_NODES];
  CommunicatorTraits::mpi_request_type recvPsi[MAX_NUM_SHARED_NODES];
  CommunicatorTraits::mpi_request_type sendPsi[MAX_NUM_SHARED_NODES];
  CommunicatorTraits::mpi_status_type statusPos[MAX_NUM_SHARED_NODES];
  CommunicatorTraits::mpi_status_type statusPsi[MAX_NUM_SHARED_NODES];
  PosType pos(P[iat]);
  //send the current position and start recv
  for(int p=0; p<NumRemoteNodes; p++)
  {
    int target=RemoteNodes[p];
    MPI_Isend(pos.data(), OHMMS_DIM, MPI_DOUBLE, target, PosTag, myComm->getMPI(), &(sendPos[p]));
    MPI_Irecv(RecvBuffer[target]->data(), PSI_DIM*OrbitalCount[target], target, PsiTag, myComm->getMPI(),&(recvPsi[p]));
  }
  ValueVector_t psiL(OrbitalCount[myNodeID]);
  GradVector_t dpsiL(OrbitalCount[myNodeID]);
  ValueVector_t p2siL(OrbitalCount[myNodeID]);
  //do the local calculation
  Phi->evaluate(pos,psiL,dpsiL,d2psiL);
  copy(psiL.begin(),psiL.end(),psi.begin()+OrbitalOffset[myNoodeID]);
  copy(dpsiL.begin(),dpsiL.end(),dpsi.begin()+OrbitalOffset[myNoodeID]);
  copy(d2psiL.begin(),d2psiL.end(),d2psi.begin()+OrbitalOffset[myNoodeID]);
  int err=MPI_Waitall(NumRemoteNodes, sendPos,statusPos);
  int ngd=OrbitalCount[myNodeID]*OHMMS_DIM;
  //do the calculation with the positions recv and send back the orbitals
  for(int p=0; p<NumRemoteNodes; p++)
  {
    int target=RemoteNodes[p];
    //matching recv with MPI_Isend
    MPI_recv(Rnow[p].data(), OHMMS_DIM, MPI_DOUBLE, target, PosTag, myNodeID,myComm->getMPI(),&(statusPos[p]));
    //do the local calculation
    Phi->evaluate(Rnow[p],psiL,dpsiL,d2psiL);
    //pack the message
    SendBuffer[p]->rewind();
    for(int i=0; i<OrbitalCount[target]; i++)
    {
      SendBuffer[p]->put(psiL[i]);
      SendBuffer[p]->put(dpsiL[i].begin(),dpsiL[i].end());
      SendBuffer[p]->put(d2dpsiL[i]);
    }
    //send the wavefunction matching ready-send with MPI_Irecv
    MPI_Irsend(SendBuffer[p]->data(), PSI_DIM*OrbitalCount[target], MPI_DOUBLE, target, PsiTag, myComm->getMPI(), &(sendPsi[p]));
  }
  int nm=NumRemoteNodes,count;
  while(nm)
  {
    int err=MPI_Testsome(NumRemoteNodes,recvPsi,&count,statusPsi);
    for(int m=0; m<count; m++)
    {
      int source=statusPsi[m].MPI_SOURCE;
      RecvBuffer[source]->rewind();
      for(int t=OrbitalOffset[source]; t<OrbitalOffset[source+1]; t++)
      {
        RecvBuffer[source]->get(psi[t]);
        RecvBuffer[source]->get(dpsi[t].begin(),dpsi[t].end());
        RecvBuffer[source]->get(d2psi[t]);
      }
    }
    nm-=count;
  }
  err=MPI_Waitall(NumRemoteNodes, sendPsi,statusPsi);
}

void DistributedSPOSet::evaluate_notranspose(const ParticleSet& P, int first, int last,
    ValueMatrix_t& logdet, GradMatrix_t& dlogdet, ValueMatrix_t& d2logdet)
{
  CommunicatorTraits::mpi_request_type sendPos[MAX_NUM_SHARED_NODES];
  CommunicatorTraits::mpi_request_type recvPsi[MAX_NUM_SHARED_NODES];
  CommunicatorTraits::mpi_request_type sendPsi[MAX_NUM_SHARED_NODES];
  CommunicatorTraits::mpi_status_type statusPos[MAX_NUM_SHARED_NODES];
  CommunicatorTraits::mpi_status_type statusPsi[MAX_NUM_SHARED_NODES];
  int nat=last-first;
  std::vector<PosType> pos(nat);
  for(int iat=first,i=0; iat<last; iat++,i++)
    pos[i]=P.R[iat];
  //send the current position and start recv
  for(int p=0; p<NumRemoteNodes; p++)
  {
    int target=RemoteNodes[p];
    MPI_Isend(&(pos[0][0]), nat*OHMMS_DIM, MPI_DOUBLE, target, PosTag, myComm->getMPI(), &(sendPos[p]));
    MPI_Irecv(RecvBuffer[target]->data(), nat*PSI_DIM*OrbitalCount[target], target, PsiTag, myComm->getMPI(),&(recvPsi[p]));
  }
  ValueVector_t psiL(OrbitalCount[myNodeID]);
  GradVector_t dpsiL(OrbitalCount[myNodeID]);
  ValueVector_t p2siL(OrbitalCount[myNodeID]);
  for(int i=0; i<nat; i++)
  {
    //do the local calculation
    Phi->evaluate_notranspose(pos[i],psiL,dpsiL,d2psiL);
    //use copy
    for(int jc=0,j=OrbitalOffset[myNoodeID]; jc<OrbitalCount[myNodeID]; jc++,j++)
    {
      //logdet(j,i)=psiL[jc];
      logdet(i,j)=psiL[jc];
      dlogdet(i,j)=dpsiL[jc];
      d2logdet(i,j)=d2psiL[jc];
    }
  }
  int err=MPI_Waitall(NumRemoteNodes, sendPos,statusPos);
  //do the calculation with the positions recv and send back the orbitals
  for(int p=0; p<NumRemoteNodes; p++)
  {
    int target=RemoteNodes[p];
    //matching recv with MPI_Isend
    MPI_recv(Rnow[0].data(), nat*OHMMS_DIM, MPI_DOUBLE, target, PosTag, myComm->getMPI(),&(statusPos[p]));
    //evaluate for the target node
    SendBuffer[p]->rewind();
    for(int i=0; i<nat; i++)
    {
      Phi->evaluate_notranspose(Rnow[i],psiL,dpsiL,d2psiL);
      for(int j=0; j<OrbitalCount[target]; j++)
      {
        SendBuffer[p]->put(psiL[j]);
        SendBuffer[p]->put(dpsiL[j].begin(),dpsiL[j].end());
        SendBuffer[p]->put(d2dpsiL[j]);
      }
    }
    //send the wavefunction matching ready-send with MPI_Irecv
    MPI_Irsend(SendBuffer[p]->data(), nat*PSI_DIM*OrbitalCount[target], MPI_DOUBLE, target, PsiTag, myComm->getMPI(), &(sendPsi[p]));
  }
  int nm=NumRemoteNodes,count;
  while(nm)
  {
    int err=MPI_Testsome(NumRemoteNodes,recvPsi,&count,statusPsi);
    for(int m=0; m<count; m++)
    {
      int source=statusPsi[m].MPI_SOURCE;
      RecvBuffer[source]->rewind();
      for(int i=0; i<nat; i++)
      {
        for(int t=OrbitalOffset[source]; t<OrbitalOffset[source+1]; t++)
        {
          RecvBuffer[source]->get(logdet(i,t));
          //RecvBuffer[source]->get(logdet(t,i));
          RecvBuffer[source]->get(dlogdet(i,t).begin(),dlogdet(i,t).end());
          RecvBuffer[source]->get(d2logdet(i,t));
        }
      }
    }
    nm-=count;
  }
  err=MPI_Waitall(NumRemoteNodes, sendPsi,statusPsi);
}
}
