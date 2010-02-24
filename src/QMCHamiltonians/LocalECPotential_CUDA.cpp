//////////////////////////////////////////////////////////////////
// (c) Copyright 2010-  by Ken Esler and Jeongnim Kim
//////////////////////////////////////////////////////////////////
//   National Center for Supercomputing Applications &
//   Materials Computation Center
//   University of Illinois, Urbana-Champaign
//   Urbana, IL 61801
//   e-mail: esler@uiuc.edu
//
// Supported by 
//   National Center for Supercomputing Applications, UIUC
//   Materials Computation Center, UIUC
//////////////////////////////////////////////////////////////////

#include "QMCHamiltonians/LocalECPotential_CUDA.h"
#include "Particle/MCWalkerConfiguration.h"

namespace qmcplusplus {
  
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
    for (int cgroup=0; cgroup<NumIonSpecies; cgroup++) {
      IonFirst.push_back(index);
      for (int i=0; i<NumIons; i++) {
	if (ions.GroupID[i] == cgroup) {
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
    RadialPotentialType* savefunc = Vspec[groupID];
    LocalECPotential::add(groupID, ppot, z);
    RadialPotentialType* rfunc = Vspec[groupID];
    if (rfunc != savefunc) {
      // Setup CUDA spline
      SRSplines[groupID] = new TextureSpline();
      int np = rfunc->size();
      vector<RealType> scaledData(np);
      for (int ir=0; ir<np; ir++)
	scaledData[ir] = -z * (*rfunc)(ir);
      SRSplines[groupID]->set
	(&scaledData[0], np, rfunc->grid().rmin(), rfunc->grid().rmax());
    }
  }
  

  void 
  LocalECPotential_CUDA::addEnergy(MCWalkerConfiguration &W, 
			       vector<RealType> &LocalEnergy)
  {
    vector<Walker_t*> &walkers = W.WalkerList;
    
    int nw = walkers.size();
    int N = NumElecs;
    if (SumGPU.size() < nw) {
      SumGPU.resize(nw);
      SumHost.resize(nw);
    }
    for (int iw=0; iw<nw; iw++) 
      SumHost[iw] = 0.0;
    SumGPU = SumHost;

    // First, do short-range part
    vector<double> esum(nw, 0.0);
    for (int sp=0; sp<NumIonSpecies; sp++) {
      if (SRSplines[sp]) {
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
    for (int iw=0; iw<walkers.size(); iw++) {
      // fprintf (stderr, "Energy = %18.6f\n", SumHost[iw]);
      walkers[iw]->getPropertyBase()[NUMPROPERTIES+myIndex] = 
	esum[iw];
      LocalEnergy[iw] += esum[iw];
    }
  }


  QMCHamiltonianBase* 
  LocalECPotential_CUDA::makeClone(ParticleSet& qp, TrialWaveFunction& psi)
  {
    LocalECPotential_CUDA* myclone =
      new LocalECPotential_CUDA(IonRef, qp);
    for(int ig=0; ig<PPset.size(); ++ig) {
      if(PPset[ig]) {
        RadialPotentialType* ppot=PPset[ig]->makeClone();
        myclone->add(ig,ppot,gZeff[ig]);
      }
    }
    return myclone;
  }

}
