//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2021 QMCPACK developers.
//
// File developed by:   Raymond Clay III, rclay@sandia.gov, Sandia National Laboratories
//
// File created by:   Raymond Clay III, rclay@sandia.gov, Sandia National Laboratories
//////////////////////////////////////////////////////////////////////////////////////

#include "QMCWaveFunctions/TWFPrototype.h"
#include <iostream>
namespace qmcplusplus
{
TWFPrototype::TWFPrototype()
{
  std::cout<<"TWFPrototype initialized\n"; 
}

TWFPrototype::IndexType TWFPrototype::get_group_index(const IndexType gid)
{
  IndexType ndets=groups.size();
  for(IndexType i=0; i<ndets; i++)
    if(gid==groups[i]) return i;

  return IndexType(-1);
}
void TWFPrototype::add_determinant(const ParticleSet& P, const IndexType gid, SPOSet* spo)
{
  app_log()<<" Gonna add a determinant with "<<gid<<" "<<spo<<std::endl;
  if( std::find(groups.begin(), groups.end(), gid) == groups.end())
  {
    app_log()<<"adding the determinant\n";
    groups.push_back(gid);
    spos.push_back(spo);
    IndexType first = P.first(gid);
    IndexType last = P.last(gid);
    IndexType norbs = spo->getOrbitalSetSize();
    num_orbs.push_back(norbs);
    num_ptcls.push_back(last-first);
  } 
}

void TWFPrototype::get_M(const ParticleSet& P, std::vector<ValueMatrix_t>& mvec)
{
  IndexType ndets=mvec.size();
  IndexType norbs=0;
  IndexType nptcls=0;
  IndexType gid=0;
  IndexType first=0;
  IndexType last=0;
  for(IndexType i=0; i<ndets; i++)
  { 
    gid=groups[i];
    first=P.first(i);
    last=P.last(i);
    nptcls=last-first;
    norbs=spos[i]->getOrbitalSetSize();
    mvec[i]=0;
    GradMatrix_t  tmpgmat;
    ValueMatrix_t tmplmat;
    tmpgmat.resize(nptcls,norbs);
    tmplmat.resize(nptcls,norbs);
    spos[i]->evaluate_notranspose(P,first,last,mvec[i],tmpgmat,tmplmat);
   
  }
}

void TWFPrototype::get_egrad_elapl_M(const ParticleSet& P, std::vector<ValueMatrix_t>& mvec, std::vector<GradMatrix_t>& gmat, std::vector<ValueMatrix_t>& lmat)
{
  IndexType ndets=mvec.size();
  IndexType norbs=0;
  IndexType nptcls=0;
  IndexType gid=0;
  IndexType first=0;
  IndexType last=0;
  for(IndexType i=0; i<ndets; i++)
  { 
    gid=groups[i];
    first=P.first(i);
    last=P.last(i);
    nptcls=last-first;
    norbs=spos[i]->getOrbitalSetSize();
    mvec[i]=0;
    gmat[i]=0;
    lmat[i]=0;
    spos[i]->evaluate_notranspose(P,first,last,mvec[i],gmat[i],lmat[i]);
  }
}

TWFPrototype::IndexType TWFPrototype::get_M_row(const ParticleSet& P, IndexType iel, ValueVector_t& val)
{
  IndexType gid = P.getGroupID(iel);
  IndexType first= P.first(gid);
  IndexType thisIndex= iel-first;
  assert(thisIndex >= 0 && thisIndex < P.getTotalNum());
 
  IndexType detIndex=-1;
  for(IndexType i=0; i<groups.size(); i++)
  {
    if(gid==groups[i]) detIndex=i;
  }
  assert(detIndex != -1);
  GradVector_t tempg;
  ValueVector_t templ;

  IndexType norbs=spos[detIndex]->getOrbitalSetSize();
  
  tempg.resize(norbs);
  templ.resize(norbs);

  spos[detIndex]->evaluateVGL(P,iel,val,tempg,templ);

  return detIndex;
}


TWFPrototype::RealType TWFPrototype::evaluateLog(ParticleSet& P)
{
  return 0;
}

}
