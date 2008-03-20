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
#ifndef QMCPLUSPLUS_DIFFERENTIAL_TWOBODYJASTROW_H
#define QMCPLUSPLUS_DIFFERENTIAL_TWOBODYJASTROW_H
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
  class DiffTwoBodyJastrowOrbital: public DiffOrbitalBase {

    ///number of variables this object handles
    int NumVars;
    ///number of target particles
    int NumPtcls;
    ///read-only distance table
    const DistanceTableData* d_table;
    ///container for the Jastrow functions  for all tghe pairs
    vector<FT*> F;

    vector<pair<int,int> > OffSet;
    Vector<RealType> dLogPsi;
    vector<GradVectorType*> gradLogPsi;
    vector<ValueVectorType*> lapLogPsi;
    std::map<std::string,FT*> J2Unique;

  public:

    ///constructor
    DiffTwoBodyJastrowOrbital(ParticleSet& p):NumVars(0) 
    {
      NumPtcls=p.getTotalNum();
      d_table=DistanceTable::add(p);

      //reserve space for the pairs
      int ng=p.groups();
      F.reserve(ng*ng);
      OffSet.reserve(ng*ng);
    }

    ~DiffTwoBodyJastrowOrbital()
    { 
      delete_iter(gradLogPsi.begin(),gradLogPsi.end());
      delete_iter(lapLogPsi.begin(),lapLogPsi.end());
    }

    ///insert a named functor
    void insert(const string& aname, FT* j) 
    {
      J2Unique[aname]=j;
      NumVars+=j->getNumOfVariables();
    }

    void insert(std::map<std::string,FT*>& j2unique) 
    {
      J2Unique.insert(j2unique.begin(),j2unique.end());
      typename std::map<std::string,FT*>::iterator it(J2Unique.begin()),it_end(J2Unique.end());
      while(it != it_end) 
      {
        NumVars+=(*it).second->size();
        ++it;
      }
    }

    ///reset the value of all the unique Two-Body Jastrow functions
    void resetParameters(OptimizableSetType& optVariables) 
    { 
      typename std::map<std::string,FT*>::iterator it(J2Unique.begin()),it_end(J2Unique.end());
      while(it != it_end) 
      {
        (*it).second->resetParameters(optVariables); 
        ++it;
      }
    }

    /** Add a functor for a pair
     */
    void addFunc(FT* afunc)
    { 
      F.push_back(afunc); 
      OffSet.push_back(pair<int,int>(afunc->FirstIndex,afunc->LastIndex));
    }

    ///reset the distance table
    void resetTargetParticleSet(ParticleSet& P) 
    {
      d_table = DistanceTable::add(P);
    }

    void initialize()
    {
      int n=LastIndex-FirstIndex;

      cout << "Index[" << FirstIndex << " " << LastIndex << endl;

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

      //correct the offset 
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

      for(int i=0; i<d_table->size(SourceIndex); ++i) {
	for(int nn=d_table->M[i]; nn<d_table->M[i+1]; ++nn) {
	  int j = d_table->J[nn];
          int ptype=d_table->PairID[nn];

          //functor evalaute value, dudr, d2udr2 at r for all the variables
          F[ptype]->evaluate(d_table->r(nn),derivs);

          RealType rinv(d_table->rinv(nn));
          PosType dr(d_table->dr(nn));

          for(int p=OffSet[ptype].first, ip=0; p<OffSet[ptype].second; ++p,++ip)
          {
            RealType dudr(rinv*derivs[ip][1]);
            RealType lap(derivs[ip][2]+2.0*dudr);
            PosType gr(dudr*dr);

            dLogPsi[p]-=derivs[ip][0];
            (*gradLogPsi[p])[i] += gr;
            (*gradLogPsi[p])[j] -= gr;
            (*lapLogPsi[p])[i] -=lap;
            (*lapLogPsi[p])[j] -=lap;
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
 * $Id: TwoBodyJastrowOrbital.h 1761 2007-02-17 23:11:59Z jnkim $ 
 ***************************************************************************/

