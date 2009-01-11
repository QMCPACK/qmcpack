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
namespace qmcplusplus 
{

  DensityEstimator::DensityEstimator(ParticleSet& elns)
  {
    for(int i=0; i<OHMMS_DIM; ++i)
      Bounds[i]=1.0/elns.Lattice.Length[i];
  }

  void DensityEstimator::resetTargetParticleSet(ParticleSet& P)
  {
  }

  DensityEstimator::Return_t DensityEstimator::evaluate(ParticleSet& P)
  {
    density=0.0;
    Return_t dmy0,dmy1,dmy2;
    for(int iat=0; iat<P.getTotalNum(); ++iat) 
    {
      PosType ru=P.Lattice.toUnit(P.R[iat]);
      int i=static_cast<int>(DeltaInv[0]*(ru[0]-std::floor(ru[0])));
      int j=static_cast<int>(DeltaInv[1]*(ru[1]-std::floor(ru[1])));
      int k=static_cast<int>(DeltaInv[2]*(ru[2]-std::floor(ru[2])));
      //int i=static_cast<int>(DeltaInv[0]*modf(ru[0],&dmy0));
      //int j=static_cast<int>(DeltaInv[1]*modf(ru[1],&dmy1));
      //int k=static_cast<int>(DeltaInv[2]*modf(ru[2],&dmy2));
      //density(static_cast<int>(Random()*DeltaInv[0])
      //    , static_cast<int>(Random()*DeltaInv[1])
      //    , static_cast<int>(Random()*DeltaInv[2])) += 1.0;
      //WRONG!!!
      density(i,j,k)+=1.0;
    }
    //do the rest
    return 0.0;
  }


  void DensityEstimator::addObservables(PropertySetType& plist)
  {
    myIndex=plist.size();
    for(int i=0; i<density.size(0); ++i)
      for(int j=0; j<density.size(1); ++j)
        for(int k=0; k<density.size(2); ++k)
        {
          ostringstream h;
          h << myName << "_" << i << "_" << j << "_" << k;
          int dum=plist.add(h.str());
        }
  }

  void DensityEstimator::registerObservables(vector<observable_helper*>& h5list
      , hid_t gid) const
  {
    int loc=h5list.size();
    //determine the shape
    vector<int> ndim(OHMMS_DIM);
    for(int i=0; i<OHMMS_DIM; ++i) ndim[i]=density.size(i);

    observable_helper* h5o=new observable_helper(myName);
    h5o->set_dimensions(ndim,myIndex);
    h5o->open(gid);
    h5list.push_back(h5o);
  }

  void DensityEstimator::setObservables(PropertySetType& plist)
  {
    std::copy(density.first_address(),density.last_address(),plist.begin()+myIndex);
  }

  void DensityEstimator::setParticlePropertyList(PropertySetType& plist
      , int offset)
  {
    std::copy(density.first_address(),density.last_address(),plist.begin()+myIndex+offset);
  }

  bool DensityEstimator::put(xmlNodePtr cur)
  {
    Delta=0.1;
    resize();
  }

  bool DensityEstimator::get(std::ostream& os) const
  {
    os << myName << " bin =" << Delta << " bohrs " << endl;
  }

  QMCHamiltonianBase* DensityEstimator::makeClone(ParticleSet& qp
      , TrialWaveFunction& psi)
  {
    //default constructor is sufficient
    return new DensityEstimator(*this);
  }

  void  DensityEstimator::resize()
  {
    TinyVector<int,OHMMS_DIM> ng;
    for(int i=0; i<OHMMS_DIM; ++i)
    {
      DeltaInv[i]=1.0/Delta[i];
      ng[i]=static_cast<int>(DeltaInv[i]);
      if(ng[i]<2) 
      {
        APP_ABORT("DensityEstimator::resize invalid bin size");
      }
    }
    app_log() << " DensityEstimator bin_size= " <<ng << " delta = " << Delta << endl;
    density.resize(ng[0],ng[1],ng[2]);
  }
}

/***************************************************************************
 * $RCSfile$   $Author: jnkim $
 * $Revision: 2945 $   $Date: 2008-08-05 10:21:33 -0500 (Tue, 05 Aug 2008) $
 * $Id: ForceBase.h 2945 2008-08-05 15:21:33Z jnkim $ 
 ***************************************************************************/
