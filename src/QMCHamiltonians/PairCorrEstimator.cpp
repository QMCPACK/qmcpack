//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Ken Esler, kpesler@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Jaron T. Krogel, krogeljt@ornl.gov, Oak Ridge National Laboratory
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////
    
    
#include <QMCHamiltonians/PairCorrEstimator.h>
#include <Particle/DistanceTableData.h>
#include <OhmmsData/AttributeSet.h>
#include <Utilities/SimpleParser.h>
#include <set>

namespace qmcplusplus
{

PairCorrEstimator::PairCorrEstimator(ParticleSet& elns, std::string& sources)
  :Dmax(10.), Delta(0.5), num_species(2)
{
  UpdateMode.set(COLLECTABLE,1);
  num_species=elns.groups();
  n_vec.resize(num_species,0);
  for(int i(0); i<num_species; i++)
    n_vec[i]=1.0/(elns.last(i)-elns.first(i));
  N_e=elns.getTotalNum();
  //use the simulation cell radius if any direction is periodic
  if(elns.Lattice.SuperCellEnum)
  {
    Dmax=elns.Lattice.WignerSeitzRadius;
    Volume=elns.Lattice.Volume;
  }
  else
    Volume=1.0; //open bcs
  NumBins=static_cast<int>(Dmax/Delta);
  Delta=Dmax/static_cast<RealType>(NumBins);
  DeltaInv=1.0/Delta;
  //ostringstream h;
  //h<<"gofr_" << elns.getName();
  //gof_r_prefix.push_back(h.str());
  std::map<int,int> pair_map;
  int npairs=0;
  for(int i=0; i<num_species; ++i)
    for(int j=i; j<num_species; ++j, ++npairs)
    {
      std::ostringstream os;
      os << "gofr_" << elns.getName() << "_" << i << "_" << j;
      gof_r_prefix.push_back(os.str());
      pair_map[i*num_species+j]=npairs;
    }
  const DistanceTableData&  dii(*elns.DistTables[0]);
  if(dii.DTType == DT_AOS)
  {
    pair_ids.resize(dii.getTotNadj());
    for(int iat=0; iat<dii.centers(); ++iat)
    {
      int ig=elns.GroupID[iat]*num_species;
      for(int nn=dii.M[iat]; nn<dii.M[iat+1]; ++nn)
        pair_ids[nn]=pair_map[ig+elns.GroupID[dii.J[nn]]];
    }
  }
  // source-target tables
  std::vector<std::string> slist, dlist;
  for(int k=1; k<elns.DistTables.size(); ++k)
    dlist.push_back(elns.DistTables[k]->origin().getName());
  parsewords(sources.c_str(),slist);
  int ntables=elns.DistTables.size();
  std::set<int> others_sorted;
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
  copy(others_sorted.begin(),others_sorted.end(),other_ids.begin());
  int toff=gof_r_prefix.size();
  for(int k=0; k<other_ids.size(); ++k)
  {
    const DistanceTableData& t(*elns.DistTables[other_ids[k]]);
    app_log() << "  GOFR for " << t.getName() << " starts at " << toff << std::endl;
    other_offsets[k]=toff;
    const SpeciesSet& species(t.origin().getSpeciesSet());
    int ng=species.size();
    for(int i=0; i<ng; ++i)
    {
      std::ostringstream os;
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
  BufferType& collectables(P.Collectables);
  const DistanceTableData&  dii(*P.DistTables[0]);
  if(dii.DTType == DT_SOA)
  {
    for(int iat=1; iat<dii.centers(); ++iat)
    {
      const RealType* restrict dist=dii.Distances[iat];
      const int ig=P.GroupID[iat];
      for(int j=0; j<iat; ++j)
      {
        const RealType r=dist[j];
        if(r>=Dmax)
          continue;
        const int loc=static_cast<int>(DeltaInv*r);
        const int jg=P.GroupID[j];
        const int pair_id=ig*(ig+1)/2+jg;
        collectables[pair_id*NumBins+loc+myIndex] += norm_factor(pair_id,loc);
      }
    }
    for(int k=0; k<other_ids.size(); ++k)
    {
      const DistanceTableData&  d1(*P.DistTables[other_ids[k]]);
      const ParticleSet::ParticleIndex_t& gid(d1.origin().GroupID);
      int koff=other_offsets[k];
      RealType overNI=1.0/d1.centers();
      for(int iat=0; iat<d1.targets(); ++iat)
      {
        const RealType* restrict dist=d1.Distances[iat];
        for(int j=0; j<d1.centers(); ++j)
        {
          const RealType r=dist[j];
          if(r>=Dmax)
            continue;
          int toff= (gid[j]+koff)*NumBins;
          int loc=static_cast<int>(DeltaInv*r);
          collectables[toff+loc+myIndex] += norm_factor(0,loc)*overNI;
        }
      }
    }
  }
  else
  {
    for(int iat=0; iat<dii.centers(); ++iat)
    {
      for(int nn=dii.M[iat]; nn<dii.M[iat+1]; ++nn)
      {
        RealType r=dii.r(nn);
        if(r>=Dmax)
          continue;
        int loc=static_cast<int>(DeltaInv*r);
        collectables[pair_ids[nn]*NumBins+loc+myIndex] += norm_factor(pair_ids[nn]+1,loc);
      }
    }
    for(int k=0; k<other_ids.size(); ++k)
    {
      const DistanceTableData&  d1(*P.DistTables[other_ids[k]]);
      const ParticleSet::ParticleIndex_t& gid(d1.origin().GroupID);
      int koff=other_offsets[k];
      for(int iat=0; iat<d1.centers(); ++iat)
      {
        RealType overNI=1.0/d1.centers();
        int toff= (gid[iat]+koff)*NumBins;
        for(int nn=d1.M[iat]; nn<d1.M[iat+1]; ++nn)
        {
          RealType r=dii.r(nn);
          if(r>=Dmax)
            continue;
          int loc=static_cast<int>(DeltaInv*r);
          collectables[toff+loc+myIndex] += norm_factor(0,loc)*overNI;
        }
      }
    }
  }
  return 0.0;
}

void PairCorrEstimator::registerCollectables(std::vector<observable_helper*>& h5list
    , hid_t gid) const
{
  std::vector<int> onedim(1,NumBins);
  int offset=myIndex;
  for(int i=0; i<gof_r_prefix.size(); ++i)
  {
    observable_helper* h5o=new observable_helper(gof_r_prefix[i]);
    h5o->set_dimensions(onedim,offset);
    h5o->open(gid);
    h5o->addProperty(const_cast<RealType&>(Delta),"delta");
    h5o->addProperty(const_cast<RealType&>(Dmax),"cutoff");
//       h5o->addProperty(const_cast<std::vector<RealType>&>(norm_factor),"norm_factor");
//       std::string blob("norm_factor[i]=1/r_m[i]^2 for r_m[i]=(r[i]+r[i+1])/2");
//       h5o->addProperty(blob,"dictionary");
    h5list.push_back(h5o);
    offset+=NumBins;
  }
}


void PairCorrEstimator::addObservables(PropertySetType& plist, BufferType& collectables)
{
  myIndex=collectables.size();
  std::vector<RealType> g(gof_r_prefix.size()*NumBins,0);
  collectables.add(g.begin(),g.end());
  ////only while debugging
  //if(gof_r.size())
  //{
  //  myDebugIndex=plist.size();
  //  for(int i=0; i<gof_r_prefix.size(); ++i)
  //  {
  //    for(int k=0; k<gof_r.cols(); ++k)
  //    {
  //      std::ostringstream h;
  //      h << gof_r_prefix[i]<< "_" << k;
  //      int dum=plist.add(h.str());
  //    }
  //  }
  //}
}


void PairCorrEstimator::setObservables(PropertySetType& plist)
{
  //std::copy(gof_r.first_address(),gof_r.last_address(),plist.begin()+myIndex);
}

void PairCorrEstimator::setParticlePropertyList(PropertySetType& plist
    , int offset)
{
  //std::copy(gof_r.first_address(),gof_r.last_address(),plist.begin()+myDebugIndex+offset);
}

bool PairCorrEstimator::put(xmlNodePtr cur)
{
  //set resolution
  int nbins=(int)std::ceil(Dmax*DeltaInv);
  std::string debug("no");
  OhmmsAttributeSet attrib;
  attrib.add(nbins,"num_bin");
  attrib.add(Dmax,"rmax");
  attrib.add(Delta,"dr");
  attrib.add(debug,"debug");
  attrib.put(cur);
  Delta=Dmax/static_cast<RealType>(nbins);
  DeltaInv=1.0/Delta;
  NumBins=nbins;
  norm_factor.resize((num_species*num_species-num_species)/2+num_species+1,NumBins);
  RealType r=Delta*0.5;
  RealType pf=Volume*DeltaInv/(4*M_PI);
  for(int i=0; i<NumBins; ++i, r+=Delta)
  {
    RealType rm2=pf/r/r;
    norm_factor(0,i)=rm2/N_e;
    int indx(1);
    for(int m(0); m<num_species; m++)
      for(int n(m); n<num_species; n++)
      {
        norm_factor(indx,i)=(m==n?2:1)*rm2*n_vec[n]*n_vec[m];
        indx++;
      }
  }
//     for(int m(0);m<norm_factor.size1();m++)
//     {
//       std::cerr <<m<<": ";
//       for(int i=0; i<NumBins; ++i)
//         std::cerr <<" "<<norm_factor(m,i);
//       std::cerr << std::endl;
//     }
  ////resize(nbins);
  //if(debug == "yes")
  //  gof_r.resize(gof_r_prefix.size(),NumBins);

  report();

  return true;
}

  
void PairCorrEstimator::report()
{
  app_log()<<"PairCorrEstimator report"<< std::endl;
  app_log()<<"  num_species = "<< num_species << std::endl;
  app_log()<<"  Volume      = "<< Volume << std::endl;
  app_log()<<"  Dmax        = "<< Dmax << std::endl;
  app_log()<<"  NumBins     = "<< NumBins << std::endl;
  app_log()<<"  Delta       = "<< Delta << std::endl;
  app_log()<<"  DeltaInv    = "<< DeltaInv << std::endl;
  //app_log()<<"  x = "<< x << std::endl;
  app_log()<<"end PairCorrEstimator report"<< std::endl;
}

bool PairCorrEstimator::get(std::ostream& os) const
{
  os << myName << " dmax=" << Dmax << std::endl;
  return true;
}

QMCHamiltonianBase* PairCorrEstimator::makeClone(ParticleSet& qp
    , TrialWaveFunction& psi)
{
  //default constructor is sufficient
  return new PairCorrEstimator(*this);
}

void  PairCorrEstimator::resize(int nbins)
{
  NumBins=nbins;
  norm_factor.resize((num_species*num_species-num_species)/2+num_species+1,NumBins);
  RealType r=Delta*0.5;
  RealType pf=Volume*DeltaInv/(4*3.14159265);
  for(int i=0; i<NumBins; ++i, r+=Delta)
  {
    RealType rm2=pf/r/r;
    norm_factor(0,i)=rm2/N_e;
    int indx(1);
    for(int m(0); m<num_species; m++)
      for(int n(m); n<num_species; n++,indx++)
        norm_factor(indx,i)=rm2*n_vec[n]*n_vec[m];
  }
//     norm_factor.resize(nbins);
//     RealType r=Delta*0.5;
//     for(int i=0; i<norm_factor.size(); ++i, r+=Delta)
//     {
//       norm_factor[i]=1.0/r/r;
//     }
}
}

