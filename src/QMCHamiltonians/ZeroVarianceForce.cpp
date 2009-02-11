#include "QMCHamiltonians/ZeroVarianceForce.h"
#include "Particle/DistanceTable.h"
#include "Particle/DistanceTableData.h"
#include "QMCWaveFunctions/TrialWaveFunction.h"
#include "ParticleBase/ParticleAttribOps.h"

namespace qmcplusplus {

  ZeroVarianceForce::ZeroVarianceForce(ParticleSet& ions, ParticleSet& elns,
			 TrialWaveFunction &psi)
    : ForceBase(ions, elns), Ions(ions), Electrons(elns),
      Psi(psi)
  {
    for (int dim=0; dim<OHMMS_DIM; dim++) {
      grad_grad_psi[dim].resize(Nel);
      lapl_grad_psi[dim].resize(Nel);
    }
    F_ZV.resize(Nnuc);
  }

  void 
  ZeroVarianceForce::resetTargetParticleSet(ParticleSet& P) {
    int tid=P.addTable(Ions);
    if(tid != myTableIndex)  
      APP_ABORT("ZeroVarianceForce::resetTargetParticleSet found inconsistent table index");
  }
  
  void 
  ZeroVarianceForce::addObservables(QMCTraits::PropertySetType& plist) {
    QMCHamiltonianBase::addObservables(plist);
    if(FirstForceIndex<0) 
      FirstForceIndex=plist.size();
    for(int iat=0; iat<Nnuc; iat++) {
      for(int x=0; x<OHMMS_DIM; x++) {
        ostringstream obsName;
        obsName << "F_ZV" << "_" << iat << "_" << x;
	app_log() << "Adding " << obsName.str() << " to observable list.\n";
        plist.add(obsName.str());
      }
    }
  }

  void 
  ZeroVarianceForce::registerObservables(vector<observable_helper*>& h5list,
					 hid_t gid) const
  {
    QMCHamiltonianBase::registerObservables(h5list, gid);
    vector<int> ndim(2);
    ndim[0]=Nnuc;
    ndim[1]=OHMMS_DIM;
    observable_helper* h5o = new observable_helper("F_ZV");
    h5o->set_dimensions(ndim,FirstForceIndex);
    h5o->open(gid);
    h5list.push_back(h5o);
  }

  void 
  ZeroVarianceForce::setObservables(QMCTraits::PropertySetType& plist) 
  {
    QMCHamiltonianBase::setObservables(plist);
    int index = FirstForceIndex;
    for(int iat=0; iat<Nnuc; iat++) {
      for(int x=0; x<OHMMS_DIM; x++) 
        plist[index++] = F_ZV[iat][x];
    }
  }
   

  void 
  ZeroVarianceForce::setParticlePropertyList
  (QMCTraits::PropertySetType& plist, int offset)
  {
    QMCHamiltonianBase::setParticlePropertyList (plist, offset);
    int index = FirstForceIndex + offset;
    for(int iat=0; iat<Nnuc; iat++) 
      for(int x=0; x<OHMMS_DIM; x++) 
        plist[index++] = F_ZV[iat][x];
  }

  ZeroVarianceForce::Return_t 
  ZeroVarianceForce::evaluate(ParticleSet& P)
  {
    RealType lapl = Sum(P.L) + Dot(P.G, P.G);
    for (int ion=0; ion < Nnuc; ion++) {
      GradType grad = 
	Psi.evalGradSource(P, Ions, ion, grad_grad_psi, lapl_grad_psi);
      for (int dim=0; dim < OHMMS_DIM; dim++) {
// 	F_ZV[ion][dim] = 0.5*(Sum(lapl_grad_psi[dim]) 
// 			      + 2.0*Dot(grad_grad_psi[dim], P.G));
	F_ZV[ion][dim] = 0.5*(Sum(lapl_grad_psi[dim]) - lapl*grad[dim]);
      }
    }
    return Value=0.0;
  }
}
