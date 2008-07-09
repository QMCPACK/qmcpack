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
#include "QMCWaveFunctions/DiffOrbitalBase.h"
#include "Particle/DistanceTableData.h"
#include "Particle/DistanceTable.h"
#include "ParticleBase/ParticleAttribOps.h"
#include "Utilities/IteratorUtility.h"


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
    ///variables handled by this orbital
    opt_variables_type myVars;
    ///container for the Jastrow functions  for all the pairs
    vector<FT*> Fs;
    ///container for the unique Jastrow functions
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
        Funique.resize(CenterRef.getSpeciesSet().size(),0);
      }

      for(int i=0; i<Fs.size(); i++) 
        if(CenterRef.GroupID[i] == source_type) Fs[i]=afunc;
      Funique[source_type]=afunc;
    }


    ///reset the value of all the unique Two-Body Jastrow functions
    void resetParameters(opt_variables_type& active) 
    { 
      for(int i=0; i<Funique.size(); ++i)
        if(Funique[i]) Funique[i]->resetParameters(active);
    }

    void checkOutVariables(const opt_variables_type& active)
    {
      myVars.clear();
      for(int i=0; i<Funique.size(); ++i)
      {
        if(Funique[i]) 
        {
          Funique[i]->myVars.getIndex(active);
          myVars.insertFrom(Funique[i]->myVars);
        }
      }
      myVars.getIndex(active);
      NumVars=myVars.size();

      //myVars.print(cout);

      if(NumVars && dLogPsi.size()==0)
      {
        dLogPsi.resize(NumVars);
        gradLogPsi.resize(NumVars,0);
        lapLogPsi.resize(NumVars,0);
        for(int i=0; i<NumVars; ++i)
        {
          gradLogPsi[i]=new GradVectorType(NumPtcls);
          lapLogPsi[i]=new ValueVectorType(NumPtcls);
        }

        OffSet.resize(Fs.size());
        int varoffset=myVars.Index[0];
        for(int i=0; i<Fs.size(); ++i)
        {
          if(Fs[i]) 
          {
            OffSet[i].first=Fs[i]->myVars.Index.front()-varoffset;
            OffSet[i].second=Fs[i]->myVars.Index.size()+OffSet[i].first;
          }
        }
      }
    }

    ///reset the distance table
    void resetTargetParticleSet(ParticleSet& P) 
    {
      d_table = DistanceTable::add(CenterRef,P);
    }

    void evaluateDerivatives(ParticleSet& P, RealType ke0, 
        const opt_variables_type& active,
        vector<RealType>& dlogpsi,
        vector<RealType>& dhpsioverpsi)
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
	for(int nn=d_table->M[i]; nn<d_table->M[i+1]; ++nn) 
        {
          std::fill(derivs.begin(),derivs.end(),0.0);
          if(!func->evaluateDerivatives(d_table->r(nn),derivs)) continue;
	  int j = d_table->J[nn];
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

      for(int k=0; k<myVars.size(); ++k)
      {
        int kk=myVars.where(k);
        if(kk<0) continue;
        dlogpsi[kk]=dLogPsi[k];
        dhpsioverpsi[kk]=-0.5*Sum(*lapLogPsi[k])-Dot(P.G,*gradLogPsi[k]);
        //optVars.setDeriv(p,dLogPsi[ip],-0.5*Sum(*lapLogPsi[ip])-Dot(P.G,*gradLogPsi[ip]));
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

