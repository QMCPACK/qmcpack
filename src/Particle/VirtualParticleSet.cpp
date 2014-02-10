//////////////////////////////////////////////////////////////////
// (c) Copyright 2013-  by Jeongnim Kim
//////////////////////////////////////////////////////////////////
// -*- C++ -*-
#include <Particle/VirtualParticleSet.h>
#include <Particle/DistanceTableData.h>
#include <Particle/DistanceTable.h>
#include <Utilities/IteratorUtility.h>

namespace qmcplusplus
{

  VirtualParticleSet::VirtualParticleSet(ParticleSet* p, int nptcl)
    :myPtcl(p)
  {
    setName("vp");
    init_minimum(nptcl);
  }

  VirtualParticleSet::~VirtualParticleSet()
  {
    //delete_iter(DistTables.begin(),DistTables.end());
  }

  void VirtualParticleSet::init_minimum(int nptcl)
  {
    //make R, ID and GroupID available
    initBase();
    Lattice = myPtcl->Lattice;
    PrimitiveLattice = myPtcl->PrimitiveLattice;
    create(nptcl);

    ratios.resize(nptcl);

    if(myPtcl->DistTables.size())
    {
      DistTables.resize(myPtcl->DistTables.size());
      DistTables[0]=createDistanceTable(*myPtcl,*this);
      for(int i=1; i<myPtcl->DistTables.size(); ++i)
        DistTables[i]=createDistanceTable(myPtcl->DistTables[i]->origin(),*this);
      for(int i=0; i<DistTables.size(); ++i)
        DistTables[i]->ID=i;
    }
  }

  void VirtualParticleSet::reset(const ParticleSet* p)
  {
    if(p->getTotalNum() != myPtcl->getTotalNum())
      APP_ABORT("VirtualParticleSet::reset Inconsistent ParticleSet size");
    myPtcl=p;
  }

  void VirtualParticleSet::makeMoves(int iat, const ParticlePos_t& displ)
  {
    activePtcl=iat;
    activeGroup=myPtcl->GroupID[iat];
    activePos=myPtcl->R[iat];
    for(int i=0; i<R.size(); ++i)
      R[i]=activePos+displ[i];
    for (int i=0; i< DistTables.size(); i++)
      DistTables[i]->evaluate(*this);
  }

  void VirtualParticleSet::validate(int iel, int k)
  {
    for (int i=0; i< DistTables.size(); i++)
    {
      const DistanceTableData& dt(*(myPtcl->DistTables[i]));
      const DistanceTableData& dt2(*DistTables[i]);
      for(int j=0; j<dt2.centers(); ++j)
      {
        int kv=j*GlobalNum+k;
        if(abs(dt.Temp[j].r1-dt2.r(kv)) > 1e-12)
          cout << "WRONG " << j << " " << dt.Temp[j].r1 << " " << dt2.r(kv) << endl;
        SingleParticlePos_t d=dt.Temp[j].dr1-dt2.dr(kv);
        if(abs(dot(d,d))>1e-12)
          cout << "WRONG " << j << " " << dt.Temp[j].dr1 << " " << dt2.dr(kv) << endl;
      }
    }
  }

}

