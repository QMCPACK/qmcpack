//////////////////////////////////////////////////////////////////
// (c) Copyright 2008-  by Jeongnim Kim
//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////
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
#include <QMCHamiltonians/LocalMomentEstimator.h>
#include <Particle/DistanceTableData.h>
#include <OhmmsData/AttributeSet.h>
#include <Utilities/SimpleParser.h>
#include <set>

namespace qmcplusplus 
{

  LocalMomentEstimator::LocalMomentEstimator(ParticleSet& elns, ParticleSet& srcs, string& sources)
  {
    int num_species=elns.groups();
    const SpeciesSet& e_species(elns.getSpeciesSet());
    int ne=e_species.size();
    int neg=elns.getTotalNum(); 
    el_id.resize(neg);
    el_nrm.resize(ne);
//     set up electron identities
    for(int i=0; i<ne; ++i)
    {
      for(int iat=0; iat<neg; ++iat)
	if(elns.GroupID[iat]==e_species.speciesName[i]) 
	{
	  el_id[iat]=i;      
	  el_nrm[i] +=1;
	}
      el_nrm[i] = 1.0/el_nrm[i];
      app_log()<<"Setting electron "<<i<<" norm to "<<el_nrm[i]<<endl;
    }
    
    int num_srcs=srcs.groups();

    //use the simulation cell radius if any direction is periodic
    if(elns.Lattice.SuperCellEnum)
      Dmax=elns.Lattice.SimulationCellRadius;

    d_table = DistanceTable::add(srcs,els);
    const SpeciesSet& species(srcs.getSpeciesSet());
    int ng=species.size();
    nag=srcs.getTotalNum();
    ion_id.resize(nag);
    for(int i=0; i<ng; ++i)
    {
      app_log()<<<<" Adding local moment estimator for "<<species.speciesName[i]<<"_"<<elns.speciesName[j]<<endl;
      sstr nm;
      for(int j(0);j<num_species;j++)
	nm<<species.speciesName[i]<<"_"<<elns.speciesName[j];
      names.push_back(nm.str());
      for(int iat=0; iat<ng; ++iat)
	if(srcs.GroupID[iat]==species.speciesName[i]) ion_id[iat]=i;
    }
    
    lm.resize(num_srcs,num_species);
    lm=0.0;
  }

  void LocalMomentEstimator::resetTargetParticleSet(ParticleSet& P)
  {
    d_table = DistanceTable::add(srcs,P);
  }

  LocalMomentEstimator::Return_t LocalMomentEstimator::evaluate(ParticleSet& P)
  {
      for(int iat=0; iat<nag; ++iat) 
      {
	int j(0);
        for(int nn=d1.M[iat]; nn<d1.M[iat+1]; ++nn,j++)
        {
          RealType r=dii.r(nn);
          if(r>=Dmax) continue;
	  lm(ion_id[iat],el_id[j]) += el_nrm[el_id[j]];
        }
      }
      
    return 0.0;
  }

  void LocalMomentEstimator::registerCollectables(vector<observable_helper*>& h5list
      , hid_t gid) const
  {
  }


  void LocalMomentEstimator::addObservables(PropertySetType& plist, BufferType& collectables)
  {
    addObservables(plist);
  }


  bool LocalMomentEstimator::put(xmlNodePtr cur)
  {
    OhmmsAttributeSet attrib;
    attrib.add(Dmax,"rcut");
    attrib.put(cur);

    return true;
  }

  bool LocalMomentEstimator::get(std::ostream& os) const
  {
    os << myName << " rcut=" << Dmax << endl;
    return true;
  }

  QMCHamiltonianBase* LocalMomentEstimator::makeClone(ParticleSet& qp
      , TrialWaveFunction& psi)
  {
    //default constructor is sufficient
    return new LocalMomentEstimator(*this);
  }

}

/***************************************************************************
 * $RCSfile$   $Author: jnkim $
 * $Revision: 2945 $   $Date: 2008-08-05 10:21:33 -0500 (Tue, 05 Aug 2008) $
 * $Id: ForceBase.h 2945 2008-08-05 15:21:33Z jnkim $ 
 ***************************************************************************/
