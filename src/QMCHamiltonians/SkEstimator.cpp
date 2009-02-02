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
#include <QMCHamiltonians/SkEstimator.h>
#include <LongRange/StructFact.h>
#include <Utilities/IteratorUtility.h>
#include <OhmmsData/AttributeSet.h>

namespace qmcplusplus 
{

  SkEstimator::SkEstimator(ParticleSet& source)
  {
    UpdateMode.set(COLLECTABLE,1);
    NumSpecies=source.getSpeciesSet().getTotalNum();
    NumK=source.SK->KLists.numk;
    OneOverN=1.0/static_cast<RealType>(source.getTotalNum());
    Kshell=source.SK->KLists.kshell;
    MaxKshell=Kshell.size()-1;
    RhokTot.resize(NumK);

    Kmag.resize(MaxKshell);
    OneOverDnk.resize(MaxKshell);
    for(int ks=0, k=0; ks<MaxKshell; ks++)
    {
      Kmag[ks]=std::sqrt(source.SK->KLists.ksq[Kshell[ks]]);
      OneOverDnk[ks]=1.0/static_cast<RealType>(Kshell[ks+1]-Kshell[ks]);
    }
  }

  void SkEstimator::resetTargetParticleSet(ParticleSet& P)
  {
  }

  SkEstimator::Return_t SkEstimator::evaluate(ParticleSet& P)
  {
    //sum over species
    std::copy(P.SK->rhok[0],P.SK->rhok[0]+NumK,RhokTot.begin());
    for(int i=1; i<NumSpecies; ++i)
      accumulate_elements(P.SK->rhok[i],P.SK->rhok[i]+NumK,RhokTot.begin());

    Vector<ComplexType>::const_iterator iit(RhokTot.begin()),iit_end(RhokTot.end());
    for(int i=myIndex;iit != iit_end;++iit,++i)
      P.Collectables[i]+=(*iit).real()*(*iit).real()+(*iit).imag()*(*iit).imag();

    return 0.0;
  }

  void SkEstimator::addObservables(PropertySetType& plist, BufferType& collectables)
  {
    myIndex=collectables.size();
    vector<RealType> tmp(NumK);
    collectables.add(tmp.begin(),tmp.end());
  }

  void SkEstimator::registerCollectables(vector<observable_helper*>& h5desc
      , hid_t gid) const
  {
    app_log() << " SkEstimator::registerCollectables " << endl;
    vector<int> ndim(1,NumK);
    observable_helper* h5o=new observable_helper(myName);
    h5o->set_dimensions(ndim,myIndex);
    h5o->open(gid);
    h5desc.push_back(h5o);
  }

  bool SkEstimator::put(xmlNodePtr cur)
  {
    //string debug("no");
    //OhmmsAttributeSet attrib;
    //attrib.add(debug,"debug");
    //attrib.put(cur);
    //if(debug == "yes") SkInst.resize(NumK);
  }

  bool SkEstimator::get(std::ostream& os) const
  {
  }

  QMCHamiltonianBase* SkEstimator::makeClone(ParticleSet& qp
      , TrialWaveFunction& psi)
  {
    //default constructor is sufficient
    return new SkEstimator(*this);
  }
}

/***************************************************************************
 * $RCSfile$   $Author: jnkim $
 * $Revision: 2945 $   $Date: 2008-08-05 10:21:33 -0500 (Tue, 05 Aug 2008) $
 * $Id: ForceBase.h 2945 2008-08-05 15:21:33Z jnkim $ 
 ***************************************************************************/
