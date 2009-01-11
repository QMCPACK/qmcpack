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
#include <OhmmsData/AttributeSet.h>
#include <Utilities/SimpleParser.h>
#include <set>

namespace qmcplusplus 
{

  PairCorrEstimator::PairCorrEstimator(ParticleSet& elns, string& sources)
    :Dmax(10.), Delta(0.5)
  {
    int num_species=elns.groups();

    //use the simulation cell radius if any direction is periodic
    if(elns.Lattice.SuperCellEnum)
    {
      Dmax=elns.Lattice.SimulationCellRadius;
      Volume=elns.Lattice.Volume;
    }

    int num_bins=static_cast<int>(Dmax/Delta);
    Delta=Dmax/static_cast<RealType>(num_bins);
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

    // source-target tables
    vector<string> slist, dlist;
    for(int k=1;k<elns.DistTables.size(); ++k)
      dlist.push_back(elns.DistTables[k]->origin().getName());

    parsewords(sources.c_str(),slist);
    int ntables=elns.DistTables.size();
    set<int> others_sorted;
    for(int i=0; i<slist.size(); ++i)
    {
      int k=0;
      while(k<dlist.size())
      {
        if(slist[i] ==dlist[k])
        {
          others_sorted.insert(k+1);
          break;
        }
        ++k;
      }
    }

    other_ids.resize(others_sorted.size());
    other_offsets.resize(others_sorted.size());
    std::copy(others_sorted.begin(),others_sorted.end(),other_ids.begin());
    int toff=gof_r_prefix.size();
    for(int k=0; k<other_ids.size(); ++k)
    {
      const DistanceTableData& t(*elns.DistTables[other_ids[k]]);
      app_log() << "  GOFR for " << t.getName() << " starts at " << toff << endl;
      other_offsets[k]=toff;
      const SpeciesSet& species(t.origin().getSpeciesSet());
      int ng=species.size();
      for(int i=0; i<ng; ++i)
      {
        ostringstream os;
        os << "gofr_" << t.getName() << "_" << species.speciesName[i];
        gof_r_prefix.push_back(os.str());
      }
      toff += ng;
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
        int loc=static_cast<int>(DeltaInv*r);
        gof_r(pair_ids[nn],loc) += norm_factor[loc];
      }
    }

    for(int k=0; k<other_ids.size(); ++k)
    {
      const DistanceTableData&  d1(*P.DistTables[other_ids[k]]);
      const ParticleSet::ParticleIndex_t& gid(d1.origin().GroupID);
      int toff=other_offsets[k];
      for(int iat=0; iat<d1.centers(); ++iat) 
      {
        RealType* gofr_ptr=gof_r[gid[iat]+toff];
        for(int nn=d1.M[iat]; nn<d1.M[iat+1]; ++nn)
        {
          RealType r=dii.r(nn);
          if(r>=Dmax) continue;
          int loc=static_cast<int>(DeltaInv*r);
          gofr_ptr[loc]+=norm_factor[loc];
        }
      }
    }

    //do the rest
    return 0.0;
  }

  void PairCorrEstimator::registerObservables(vector<observable_helper*>& h5list
      , hid_t gid) const
  {
    vector<int> onedim(1,gof_r.cols());
    int offset=myIndex;

    for(int i=0; i<gof_r_prefix.size(); ++i)
    {
      observable_helper* h5o=new observable_helper(gof_r_prefix[i]);
      h5o->set_dimensions(onedim,offset);
      h5o->open(gid);

      h5o->addProperty(const_cast<RealType&>(Delta),"delta");
      h5o->addProperty(const_cast<RealType&>(Dmax),"cutoff");
      h5o->addProperty(const_cast<vector<RealType>&>(norm_factor),"norm_factor");

      std::string blob("norm_factor[i]=1/r_m[i]^2 for r_m[i]=(r[i]+r[i+1])/2");
      h5o->addProperty(blob,"dictionary");

      h5list.push_back(h5o);
      offset+=gof_r.cols();
    }
  }

  void PairCorrEstimator::addObservables(PropertySetType& plist)
  {
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
    std::copy(gof_r.first_address(),gof_r.last_address(),plist.begin()+myIndex+offset);
  }

  bool PairCorrEstimator::put(xmlNodePtr cur)
  {
    //set resolution 
    int nbins=Dmax*DeltaInv;
    OhmmsAttributeSet attrib;
    attrib.add(nbins,"num_bin");
    attrib.add(Dmax,"rmax");
    attrib.add(Delta,"dr");
    attrib.put(cur);

    Delta=Dmax/static_cast<RealType>(nbins);
    DeltaInv=1.0/Delta;
    resize(nbins);
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

  void  PairCorrEstimator::resize(int nbins)
  {
    gof_r.resize(gof_r_prefix.size(),nbins);
    norm_factor.resize(gof_r.cols());
    RealType r=Delta*0.5;
    for(int i=0; i<norm_factor.size(); ++i, r+=Delta)
    {
      norm_factor[i]=1.0/r/r;
    }
  }
}

/***************************************************************************
 * $RCSfile$   $Author: jnkim $
 * $Revision: 2945 $   $Date: 2008-08-05 10:21:33 -0500 (Tue, 05 Aug 2008) $
 * $Id: ForceBase.h 2945 2008-08-05 15:21:33Z jnkim $ 
 ***************************************************************************/
