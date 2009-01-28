#include "QMCHamiltonians/PulayForce.h"
#include "Particle/DistanceTable.h"
#include "Particle/DistanceTableData.h"

namespace qmcplusplus {

  PulayForce::PulayForce(ParticleSet& ions, ParticleSet& elns)
    : ForceBase(ions, elns), Ions(ions), Electrons(elns)
    {
      GradLogPsi.resize(Nnuc, PosType());
      EGradLogPsi.resize(Nnuc, PosType());
      WarpNorm.resize(Nelec, 0.0);
    }

  void 
  PulayForce::addObservables(QMCTraits::PropertySetType& plist) {
    QMCHamiltonianBasee::add
    if(FirstForceIndex<0) 
      FirstForceIndex=plist.size();
    for(int iat=0; iat<Nnuc; iat++) {
      for(int x=0; x<OHMMS_DIM; x++) {
        ostringstream obsName1, obsName2;
        obsName1 << "grad_log_psi" << "_" << iat << "_" << x;
        plist.add(obsName1.str());
        obsName2 << "Egrad_log_psi" << "_" << iat << "_" << x;
        plist.add(obsName2.str());
      }
    }
  }

  void 
  PulayForce::registerObservables(vector<observable_helper*>& h5list,
				  hid_t gid) const
  {
    QMCHamiltonianBase::registerObservables(h5list, gid);
    vector<int> ndim(2);
    ndim[0]=Nnuc;
    ndim[1]=OHMMS_DIM;
    observable_helper* h5o1 = new observable_helper("grad_log_psi");
    h5o1->set_dimensions(ndim,FirstForceIndex);
    h5o1->open(gid);
    h5list.push_back(h5o1);

    observable_helper* h5o2 = new observable_helper("Egrad_log_psi");
    h5o2->set_dimensions(ndim,FirstForceIndex+Nnuc*OHMMS_DIM);
    h5o2->open(gid);
    h5list.push_back(h5o2);
  }

  void 
  PulayForce::setObservables(QMCTraits::PropertySetType& plist) 
  {
    QMCHamiltonianBase::setObservables(plist);
    int index = FirstForceIndex;
    for(int iat=0; iat<Nnuc; iat++) {
      for(int x=0; x<OHMMS_DIM; x++) {
        plist[index++] = GradLogPsi[iat][x];
	plist[index++] = EGradLogPsi[iat][x];
      }
    }
  }
   

  void 
  PulayForce::setParticleSet(QMCTraits::PropertySetType& plist, int offset)
  {
    QMCHamiltonianBase::setParticleSet (plist, offset);
    int index = FirstForceIndex + offset;
    for(int iat=0; iat<Nnuc; iat++) {
      for(int x=0; x<OHMMS_DIM; x++) {
        plist[index++] = GradLogPsi[iat][x];
        plist[index++] = EGradLogPsi[iat][x];
      }
    }
  }

  void
  PulayForce::evaluate (ParticleSet &P) 
  {
    const DistanceTableData* d_ab=P.DistTables[myTableIndex];

    // Compute normalization of the warp tranform for each electron
    for (int elec=0; elec<Nel; elec++) 
      WarpNorm[elec] = 0.0;
    for (int ion=0; ion<Nnuc; ion++)
      for(int nn=d_ab.M[ion], elec=0; nn<d_ab.M[ion+1]; ++nn,++elec) 
	WarpNorm[elec] += WarpFunction(d_ab.r(nn));
    for (int elec=0; elec<Nel; elec++) 
      WarpNorm[elec] = 1.0/WarpNorm[elec];
    for (int ion=0; ion < Nnuc; ion++) {
      GradLogPsi[ion] = EGradLogPsi[ion] = PosType();
      for(int nn=d_ab.M[ion], elec=0; nn<d_ab.M[ion+1]; ++nn,++elec) 
	GradLogPsi[ion] += WarpNorm[elec] * WarpFunction(dab[nn].r(nn))
	  * P.G[elec];
      EGradLogPsi[ion] = P.Properties[LOCALENERGY] * GradLogPsi[ion];
    }

  }
