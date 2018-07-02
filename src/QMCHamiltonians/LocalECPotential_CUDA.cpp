//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Ken Esler, kpesler@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Ken Esler, kpesler@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////
    
    

#include "QMCHamiltonians/LocalECPotential_CUDA.h"
#include "Particle/MCWalkerConfiguration.h"

namespace qmcplusplus
{

LocalECPotential_CUDA::LocalECPotential_CUDA
(ParticleSet& ions, ParticleSet& elns) :
  LocalECPotential(ions,elns),
  ElecRef(elns), IonRef(ions),
  SumGPU("LocalECPotential::SumGPU"),
  IGPU("LocalECPotential::IGPU")
{
  SpeciesSet &sSet = ions.getSpeciesSet();
  NumIonSpecies = sSet.getTotalNum();
  NumIons  = ions.getTotalNum();
  NumElecs = elns.getTotalNum();
#ifdef QMC_CUDA
  // Copy center positions to GPU, sorting by GroupID
  gpu::host_vector<CUDA_PRECISION> I_host(OHMMS_DIM*NumIons);
  int index=0;
  for (int cgroup=0; cgroup<NumIonSpecies; cgroup++)
  {
    IonFirst.push_back(index);
    for (int i=0; i<NumIons; i++)
    {
      if (ions.GroupID[i] == cgroup)
      {
        for (int dim=0; dim<OHMMS_DIM; dim++)
          I_host[OHMMS_DIM*index+dim] = ions.R[i][dim];
        SortedIons.push_back(ions.R[i]);
        index++;
      }
    }
    IonLast.push_back(index-1);
  }
  IGPU = I_host;
  SRSplines.resize(NumIonSpecies,0);
#endif
}

void
LocalECPotential_CUDA::add(int groupID, RadialPotentialType* ppot, RealType z)
{
  RadialPotentialType* savefunc = PPset[groupID];
  LocalECPotential::add(groupID, ppot, z);
  RadialPotentialType* rfunc = PPset[groupID];
  if (rfunc != savefunc)
  {
    // Setup CUDA spline
    SRSplines[groupID] = new TextureSpline();
    //      int np = 10001;
    //       RealType rmax(20.0);
    //       char fname[100];
    //       snprintf (fname, 100, "local_ecp_%d.dat", groupID);
    //       FILE *fout = fopen (fname, "w");
    int np = rfunc->size();
    std::vector<RealType> scaledData(np);
    for (int ir=0; ir<np; ir++)
    {
      // double r = ((RealType)ir / (RealType)(np-1)) * rmax ;
      // 	scaledData[ir] = -z* rfunc->splint(r);
      //	fprintf (stderr, "V(%1.5f) = %1.8f\n", r, scaledData[ir]);
      //	fprintf (fout, "%16.10f %18.10e\n", r,
      // 	scaledData[ir]);
      scaledData[ir] = -z * (*rfunc)(ir);
    }
    // fclose(fout);
    //SRSplines[groupID]->set(&(scaledData[0]), np, 0.0, rmax);
    SRSplines[groupID]->set(&(scaledData[0]), rfunc->size(),
                            rfunc->grid().rmin(), rfunc->grid().rmax());
  }
}


void
LocalECPotential_CUDA::addEnergy(MCWalkerConfiguration &W,
                                 std::vector<RealType> &LocalEnergy)
{
  std::vector<Walker_t*> &walkers = W.WalkerList;
  int nw = walkers.size();
  int N = NumElecs;
  if (SumGPU.size() < nw)
  {
    SumGPU.resize(nw);
    SumHost.resize(nw);
  }
  for (int iw=0; iw<nw; iw++)
    SumHost[iw] = 0.0;
  SumGPU = SumHost;
  // First, do short-range part
  std::vector<double> esum(nw, 0.0);
  for (int sp=0; sp<NumIonSpecies; sp++)
  {
    if (SRSplines[sp])
    {
      local_ecp_sum
      (W.RList_GPU.data(), N, IGPU.data(), IonFirst[sp], IonLast[sp],
       SRSplines[sp]->rMax, SRSplines[sp]->NumPoints,
       SRSplines[sp]->MyTexture, SumGPU.data(), nw);
      SumHost = SumGPU;
      for (int iw=0; iw<nw; iw++)
        esum[iw] += SumHost[iw];
    }
  }
  SumHost = SumGPU;
  for (int iw=0; iw<walkers.size(); iw++)
  {
    // fprintf (stderr, "Energy = %18.6f\n", SumHost[iw]);
    walkers[iw]->getPropertyBase()[NUMPROPERTIES+myIndex] = esum[iw];
    LocalEnergy[iw] += esum[iw];
  }
}


QMCHamiltonianBase*
LocalECPotential_CUDA::makeClone(ParticleSet& qp, TrialWaveFunction& psi)
{
  LocalECPotential_CUDA* myclone =
    new LocalECPotential_CUDA(IonRef, qp);
  for(int ig=0; ig<PPset.size(); ++ig)
  {
    if(PPset[ig])
    {
      RadialPotentialType* ppot=PPset[ig]->makeClone();
      myclone->add(ig,ppot,gZeff[ig]);
    }
  }
  return myclone;
}

}
