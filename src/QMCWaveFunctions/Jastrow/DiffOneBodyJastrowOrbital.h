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
#ifndef QMCPLUSPLUS_DIFFERENTIAL_ONEBODYJASTROW_H
#define QMCPLUSPLUS_DIFFERENTIAL_ONEBODYJASTROW_H
#include "Configuration.h"
#include  <map>
#include "QMCWaveFunctions/DiffOrbitalBase.h"
#include "Particle/DistanceTableData.h"
#include "Particle/DistanceTable.h"
#include "ParticleBase/ParticleAttribOps.h"

namespace qmcplusplus {

  /** @ingroup OrbitalComponent
   *  @brief Specialization for two-body Jastrow function using multiple functors
   */ 
  template<class FT>
  class DiffOneBodyJastrowOrbital: public DiffOrbitalBase {

    ///number of variables this object handles
    int NumVars;
    ///number of target particles
    int NumPtcls;
    ///reference to the ions
    const ParticleSet& CenterRef;
    ///read-only distance table
    const DistanceTableData* d_table;
    ///container for the Jastrow functions  for all tghe pairs
    vector<FT*> Fs;
    vector<FT*> Funique;

    vector<pair<int,int> > OffSet;
    Vector<RealType> dLogPsi;
    vector<GradVectorType*> gradLogPsi;
    vector<ValueVectorType*> lapLogPsi;

  public:

    ///constructor
    DiffOneBodyJastrowOrbital(const ParticleSet& centers, ParticleSet& els)
      :CenterRef(centers),NumVars(0) 
    {
      NumPtcls=els.getTotalNum();
      d_table=DistanceTable::add(centers,els);
    }

    ~DiffOneBodyJastrowOrbital()
    { 
      delete_iter(gradLogPsi.begin(),gradLogPsi.end());
      delete_iter(lapLogPsi.begin(),lapLogPsi.end());
    }

    /** Add a radial functor for a group
     * @param source_type group index of the center species
     * @param afunc radial functor
     */
    void addFunc(int source_type, FT* afunc) 
    {
      if(Fs.empty()) 
      {
        Fs.resize(CenterRef.getTotalNum(),0);
        OffSet.resize(CenterRef.getTotalNum());
      }

      bool foundgroup=false;
      for(int i=0; i<Fs.size(); i++) 
      {
        if(CenterRef.GroupID[i] == source_type) 
        {
          Fs[i]=afunc;
          OffSet[i].first=afunc->FirstIndex;
          OffSet[i].second=afunc->LastIndex;
          foundgroup=true;
        }
      }

      if(foundgroup)
      {
        Funique.push_back(afunc);
        NumVars+=afunc->getNumOfVariables();
      }
    }


    ///reset the value of all the unique Two-Body Jastrow functions
    void resetParameters(OptimizableSetType& optVariables) 
    { 
      for(int i=0; i<Funique.size(); ++i)
        Funique[i]->resetParameters(optVariables);
    }

    ///reset the distance table
    void resetTargetParticleSet(ParticleSet& P) 
    {
      d_table = DistanceTable::add(CenterRef,P);
    }

    void initialize()
    {
      int n=LastIndex-FirstIndex;

      if(dLogPsi.size()==0)
      {
        dLogPsi.resize(n);
        gradLogPsi.resize(n,0);
        lapLogPsi.resize(n,0);
        for(int i=0; i<n; ++i)
        {
          gradLogPsi[i]=new GradVectorType(NumPtcls);
          lapLogPsi[i]=new ValueVectorType(NumPtcls);
        }
      }

      //correct the offset , don't need to worry about null functor
      for(int p=0; p<OffSet.size(); ++p)
      {
        OffSet[p].first-=FirstIndex;
        OffSet[p].second-=FirstIndex;
      }
    }

    void evaluateDerivatives(ParticleSet& P, RealType ke0, OptimizableSetType& optVars)
    {

      dLogPsi=0.0;
      for(int p=0;p<NumVars; ++p) (*gradLogPsi[p])=0.0;
      for(int p=0;p<NumVars; ++p) (*lapLogPsi[p])=0.0;
      
      vector<TinyVector<RealType,3> > derivs(NumVars);

      for(int i=0; i<d_table->size(SourceIndex); ++i) 
      {
        FT* func=Fs[i];
        if(func == 0) continue;
        int first(OffSet[i].first);
        int last(OffSet[i].second);
	for(int nn=d_table->M[i]; nn<d_table->M[i+1]; ++nn) {
	  int j = d_table->J[nn];
          func->evaluate(d_table->r(nn),derivs);
          RealType rinv(d_table->rinv(nn));
          PosType dr(d_table->dr(nn));
          for(int p=first, ip=0; p<last; ++p,++ip)
          {
            dLogPsi[p] -= derivs[ip][0];
            RealType dudr(rinv*derivs[ip][1]);
            (*gradLogPsi[p])[j] -= dudr*dr;
            (*lapLogPsi[p])[j] -= derivs[ip][2]+2.0*dudr;
          }
	}
      }

      for(int p=FirstIndex, ip=0; p<LastIndex; ++p,++ip)
      {
        optVars.setDeriv(p,dLogPsi[ip],-0.5*Sum(*lapLogPsi[ip])-Dot(P.G,*gradLogPsi[ip]));
      }
    }
  };
}
#endif
/***************************************************************************
 * $RCSfile$   $Author: jnkim $
 * $Revision: 1761 $   $Date: 2007-02-17 17:11:59 -0600 (Sat, 17 Feb 2007) $
 * $Id: OneBodyJastrowOrbital.h 1761 2007-02-17 23:11:59Z jnkim $ 
 ***************************************************************************/

