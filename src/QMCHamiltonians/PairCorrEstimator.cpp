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
#include <QMCHamiltonians/PairCorrEstimator.h>
#include <Particle/DistanceTableData.h>

namespace qmcplusplus 
{

  PairCorrEstimator::PairCorrEstimator(ParticleSet& elns)
    :Dmax(10.), Delta(0.5)
  {
    int num_species=elns.groups();

    //use the simulation cell radius if any direction is periodic
    if(elns.Lattice.SuperCellEnum)
      Dmax=elns.Lattice.SimulationCellRadius;
    DeltaInv=1.0/Delta;

    //ostringstream h;
    //h<<"gofr_" << elns.getName();
    //gof_r_prefix.push_back(h.str());

    map<int,int> pair_map;
    int npairs=0;
    for(int i=0; i<num_species; ++i)
      for(int j=i; j<num_species; ++j, ++npairs)
      {
        ostringstream os;
        os << "gofr_" << elns.getName() << "_" << i << "_" << j;
        gof_r_prefix.push_back(os.str());
        pair_map[i*num_species+j]=npairs;
      }

    const DistanceTableData&  dii(*elns.DistTables[0]);
    pair_ids.resize(dii.getTotNadj());
    for(int iat=0; iat<dii.centers(); ++iat) 
    {
      int ig=elns.GroupID[iat]*num_species;
      for(int nn=dii.M[iat]; nn<dii.M[iat+1]; ++nn)
        pair_ids[nn]=pair_map[ig+elns.GroupID[dii.J[nn]]];
    }
  }

  void PairCorrEstimator::resetTargetParticleSet(ParticleSet& P)
  {
  }

  PairCorrEstimator::Return_t PairCorrEstimator::evaluate(ParticleSet& P)
  {
    gof_r=0.0;
    const DistanceTableData&  dii(*P.DistTables[0]);
    for(int iat=0; iat<dii.centers(); ++iat) 
    {
      for(int nn=dii.M[iat]; nn<dii.M[iat+1]; ++nn)
      {
        RealType r=dii.r(nn);
        if(r>=Dmax) continue;
        gof_r(pair_ids[nn],static_cast<int>(DeltaInv*r)) += 1.0;
      }
    }

    //do the rest
    return 0.0;
  }


  void PairCorrEstimator::addObservables(PropertySetType& plist)
  {
    //this is an anchor
    myIndex=plist.size();
    for(int i=0; i<gof_r_prefix.size(); ++i)
    {
      for(int k=0; k<gof_r.cols(); ++k)
      {
        ostringstream h;
        h << gof_r_prefix[i]<< "_" << k;
        int dum=plist.add(h.str());
      }
    }
  }

  void PairCorrEstimator::setObservables(PropertySetType& plist)
  {
    std::copy(gof_r.first_address(),gof_r.last_address(),plist.begin()+myIndex);
  }

  void PairCorrEstimator::setParticlePropertyList(PropertySetType& plist
      , int offset)
  {
    std::copy(gof_r.first_address(),gof_r.last_address(),plist.begin()+offset);
  }

  bool PairCorrEstimator::put(xmlNodePtr cur)
  {
    resize();
  }

  bool PairCorrEstimator::get(std::ostream& os) const
  {
    os << myName << " dmax=" << Dmax << endl;
  }

  QMCHamiltonianBase* PairCorrEstimator::makeClone(ParticleSet& qp
      , TrialWaveFunction& psi)
  {
    //default constructor is sufficient
    return new PairCorrEstimator(*this);
  }

  void  PairCorrEstimator::resize()
  {
    gof_r.resize(gof_r_prefix.size(),Dmax*DeltaInv+1);
  }
}

/***************************************************************************
 * $RCSfile$   $Author: jnkim $
 * $Revision: 2945 $   $Date: 2008-08-05 10:21:33 -0500 (Tue, 05 Aug 2008) $
 * $Id: ForceBase.h 2945 2008-08-05 15:21:33Z jnkim $ 
 ***************************************************************************/
