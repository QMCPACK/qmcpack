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
#include <Particle/DistanceTableData.h>

namespace qmcplusplus 
{

  DensityEstimator::DensityEstimator(ParticleSet& elns)
  {
    for(int i=0; i<OHMMS_DIM; ++i)
      Bounds[i]=elns.Lattice.Length[i];
  }

  void DensityEstimator::resetTargetParticleSet(ParticleSet& P)
  {
  }

  DensityEstimator::Return_t DensityEstimator::evaluate(ParticleSet& P)
  {
    density=0.0;
    //not correct
    for(int iat=0; iat<P.getTotalNum(); ++iat) 
    {
      //WRONG!!!
      density(static_cast<int>(DeltaInv[0]*P.R[iat][0])
          , static_cast<int>(DeltaInv[1]*P.R[iat][1])
          , static_cast<int>(DeltaInv[2]*P.R[iat][2])) += 1.0;
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
          h << "den_" << i << "_" << j << "_" << k;
          int dum=plist.add(h.str());
        }
  }

  void DensityEstimator::setObservables(PropertySetType& plist)
  {
    std::copy(density.first_address(),density.last_address(),plist.begin()+myIndex);
  }

  void DensityEstimator::setParticlePropertyList(PropertySetType& plist
      , int offset)
  {
    std::copy(density.first_address(),density.last_address(),plist.begin()+offset);
  }

  bool DensityEstimator::put(xmlNodePtr cur)
  {
    Delta=0.5;
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
      DeltaInv[i]=Delta[i]/Bounds[i];
      ng[i]=static_cast<int>(DeltaInv[i]);
      if(ng[i]) APP_ABORT("DensityEstimator::resize invalid bin size");
    }
    density.resize(ng[0],ng[1],ng[2]);
  }
}

/***************************************************************************
 * $RCSfile$   $Author: jnkim $
 * $Revision: 2945 $   $Date: 2008-08-05 10:21:33 -0500 (Tue, 05 Aug 2008) $
 * $Id: ForceBase.h 2945 2008-08-05 15:21:33Z jnkim $ 
 ***************************************************************************/
