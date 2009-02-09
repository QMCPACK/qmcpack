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
#include <QMCHamiltonians/DensityEstimator.h>
#include <OhmmsData/AttributeSet.h>
namespace qmcplusplus 
{

  DensityEstimator::DensityEstimator(ParticleSet& elns)
  {
    UpdateMode.set(COLLECTABLE,1);
    for(int i=0; i<OHMMS_DIM; ++i)
      Bounds[i]=1.0/elns.Lattice.Length[i];
  }

  void DensityEstimator::resetTargetParticleSet(ParticleSet& P)
  {
  }

  DensityEstimator::Return_t DensityEstimator::evaluate(ParticleSet& P)
  {
    for(int iat=0; iat<P.getTotalNum(); ++iat) 
    {
      PosType ru=P.Lattice.toUnit(P.R[iat]);
      int i=static_cast<int>(DeltaInv[0]*(ru[0]-std::floor(ru[0])));
      int j=static_cast<int>(DeltaInv[1]*(ru[1]-std::floor(ru[1])));
      int k=static_cast<int>(DeltaInv[2]*(ru[2]-std::floor(ru[2])));
      P.Collectables[getGridIndex(i,j,k)]+=1.0;
    }
    return 0.0;
  }

  void DensityEstimator::addObservables(PropertySetType& plist, BufferType& collectables)
  {
    //current index
    myIndex=collectables.current();
    vector<RealType> tmp(NumGrids[OHMMS_DIM]);
    collectables.add(tmp.begin(),tmp.end());
  }

  void DensityEstimator::registerCollectables(vector<observable_helper*>& h5desc
      , hid_t gid) const
  {
    int loc=h5desc.size();
    vector<int> ng(OHMMS_DIM);
    for(int i=0; i<OHMMS_DIM; ++i) ng[i]=NumGrids[i];
    observable_helper* h5o=new observable_helper(myName);
    h5o->set_dimensions(ng,myIndex);
    h5o->open(gid);
    h5desc.push_back(h5o);
  }

  void DensityEstimator::setObservables(PropertySetType& plist)
  {
    //std::copy(density.first_address(),density.last_address(),plist.begin()+myDebugIndex);
  }

  void DensityEstimator::setParticlePropertyList(PropertySetType& plist
      , int offset)
  {
    //std::copy(density.first_address(),density.last_address(),plist.begin()+myDebugIndex+offset);
  }

  /** check xml elements
   *
   * <estimator name="density" debug="no" delta="0.1"/>
   */
  bool DensityEstimator::put(xmlNodePtr cur)
  {
    Delta=0.1;

    string debug("no");
    OhmmsAttributeSet attrib;
    attrib.add(debug,"debug");
    attrib.add(Delta,"delta");
    attrib.put(cur);

    resize();
    return true;
  }

  bool DensityEstimator::get(std::ostream& os) const
  {
    os << myName << " bin =" << Delta << " bohrs " << endl;
    return true;
  }

  QMCHamiltonianBase* DensityEstimator::makeClone(ParticleSet& qp
      , TrialWaveFunction& psi)
  {
    //default constructor is sufficient
    return new DensityEstimator(*this);
  }

  void  DensityEstimator::resize()
  {
    for(int i=0; i<OHMMS_DIM; ++i)
    {
      DeltaInv[i]=1.0/Delta[i];
      NumGrids[i]=static_cast<int>(DeltaInv[i]);
      if(NumGrids[i]<2) 
      {
        APP_ABORT("DensityEstimator::resize invalid bin size");
      }
    }
    app_log() << " DensityEstimator bin_size= " <<NumGrids << " delta = " << Delta << endl;
    NumGrids[OHMMS_DIM]=NumGrids[0]*NumGrids[1]*NumGrids[2];
  }
}

/***************************************************************************
 * $RCSfile$   $Author: jnkim $
 * $Revision: 2945 $   $Date: 2008-08-05 10:21:33 -0500 (Tue, 05 Aug 2008) $
 * $Id: ForceBase.h 2945 2008-08-05 15:21:33Z jnkim $ 
 ***************************************************************************/
