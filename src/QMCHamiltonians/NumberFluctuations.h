//////////////////////////////////////////////////////////////////
// (c) Copyright 2003  by Jeongnim Kim
//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////
//   National Center for Supercomputing Applications &
//   Materials Computation Center
//   University of Illinois, Urbana-Champaign
//   Urbana, IL 61801
//   e-mail: jnkim@ncsa.uiuc.edu
//   Tel:    217-244-6319 (NCSA) 217-333-3324 (MCC)
//
// Supported by 
//   National Center for Supercomputing Applications, UIUC
//   Materials Computation Center, UIUC
//////////////////////////////////////////////////////////////////
// -*- C++ -*-
#ifndef QMCPLUSPLUS_NUMBERFLUC_H
#define QMCPLUSPLUS_NUMBERFLUC_H
#include "Particle/ParticleSet.h"
#include "Particle/WalkerSetRef.h"
#include "Particle/DistanceTableData.h"
#include "Particle/DistanceTable.h"
#include "QMCHamiltonians/QMCHamiltonianBase.h"
#include <numeric>

namespace qmcplusplus {

  struct NumberFluctuations: public QMCHamiltonianBase {

    vector<Return_t> regions;
    string shape;
    ParticleSet::ParticlePos_t* L;
    bool PBCType;
    int nreg,grps,pwrs;
    vector<Return_t>  Values, voffsets;
    
    
    NumberFluctuations(ParticleSet& P) 
    {
      L = new ParticleSet::ParticlePos_t(P.R.size());
      L->setUnit(PosUnit::LatticeUnit);
      grps=P.groups();
      PBCType=(P.Lattice.SuperCellEnum>1);
    }

    ~NumberFluctuations() { }

    void resetTargetParticleSet(ParticleSet& P)  
    {
      grps=P.groups();
      Values.resize( (grps+1)*nreg*pwrs,0 );
    }

    inline Return_t 
    evaluateSphere(ParticleSet& P) 
    { 
      (*L)=P.R;
//       P.applyBC(*L);
      P.applyMinimumImage(*L);
      Matrix<int> np(grps+1,nreg);
      np=0;
      for (int g=0; g<grps; ++g)
        for (int iat=P.first(g); iat<P.last(g); ++iat)
        { 
          Return_t z=dot((*L)[iat],(*L)[iat]);
          int a(nreg-1);
          while((z<regions[a])&&(a>=0))
          {
            np(g,a)+=1;
            np(grps,a)+=1;
            a--;
          }
        }
       
       for (int g=0; g<Values.size(); ++g) Values[g]=0;
       
       int indx(0);
       for (int g=0; g<grps+1; ++g) for (int h=0; h<nreg; h++)
       {
         Return_t vx=np(g,h)-voffsets[g*nreg+h];
         Return_t vy(1);
         for (int i=0; i<pwrs; i++,++indx)
           Values[indx]+=(vy*=vx);
       }
       
       return 0;
    }

    inline Return_t 
    evaluateBox(ParticleSet& P) 
    { 
      P.convert2CartInBox(P.R,*L);
      Matrix<int> np(grps+1,nreg);
      np=0;
      for (int g=0; g<grps; ++g)
        for (int iat=P.first(g); iat<P.last(g); ++iat)
        { 
          Return_t z=(*L)[iat][0];
          for (int x=1;x<DIM;x++) z=std::max(z,(*L)[iat][x]);
          int a(nreg-1);
          while((z<regions[a])&&(a>=0))
          {
            np(g,a)+=1;
            np(grps,a)+=1;
            a--;
          }
        }
       
       for (int g=0; g<Values.size(); ++g) Values[g]=0;
       
       int indx(0);
       for (int g=0; g<grps+1; ++g) for (int h=0; h<nreg; h++) for (int i=0; i<pwrs; i++,++indx)
         Values[indx]+=std::pow(np(g,h)-voffsets[g*nreg+h],i+1);
       
       return 0;
    }

    inline Return_t 
    evaluateHalfSpace(ParticleSet& P) 
    { 
      if(PBCType) P.convert2CartInBox(P.R,*L);
      Matrix<int> np(grps+1,nreg);
      np=0;
      for (int g=0; g<grps; ++g)
        for (int iat=P.first(g); iat<P.last(g); ++iat)
        { 
          Return_t z=(*L)[iat][DIM-1];
          int a(nreg-1);
          while((z<regions[a])&&(a>=0))
          {
            np(g,a)+=1;
            np(grps,a)+=1;
            a--;
          }
        }
       
       for (int g=0; g<Values.size(); ++g) Values[g]=0;
       
       int indx(0);
       for (int g=0; g<grps+1; ++g) for (int h=0; h<nreg; h++) for (int i=0; i<pwrs; i++,++indx)
         Values[indx]+=std::pow(np(g,h)-voffsets[g*nreg+h],i+1);
       
       return 0;
    }
    
    void addObservables(PropertySetType& plist, BufferType& collectables)
    {
      myIndex=plist.size();
      for (int g=0; g<grps+1; ++g) for (int h=0; h<nreg; h++) for (int i=0; i<pwrs; i++)
      {
        std::stringstream sstr;
        sstr << "N"<<g<<"_R"<<h<<"_P"<<i+1;
        int id=plist.add(sstr.str());
      }
    }
    
    void setObservables(PropertySetType& plist)
    {
      std::copy(Values.begin(),Values.end(),plist.begin()+myIndex);
    }

    void setParticlePropertyList(PropertySetType& plist, int offset)
    {
      std::copy(Values.begin(),Values.end(),plist.begin()+myIndex+offset);
    }
    
    inline Return_t evaluate(ParticleSet& P)
    {
      if (shape=="sphere") return evaluateSphere(P);
      else if (shape=="box") return evaluateBox(P);
      else if (shape=="half") return evaluateHalfSpace(P);
    }

    inline Return_t evaluate(ParticleSet& P, vector<NonLocalData>& Txy) {
      return evaluate(P);
    }

    inline Return_t 
    registerData(ParticleSet& P, BufferType& buffer) 
    {
    }

    inline Return_t 
    updateBuffer(ParticleSet& P, BufferType& buffer) 
    {
    }

    inline void copyFromBuffer(ParticleSet& P, BufferType& buffer)
    {
    }

    inline void copyToBuffer(ParticleSet& P, BufferType& buffer)
    {
    }

    inline Return_t 
    evaluatePbyP(ParticleSet& P, int active)
    {
      APP_ABORT("NumberFluctuations::evaluatePbyP");
      return 0.0;
    }

    /** Do nothing */
    bool put(xmlNodePtr cur) {
      OhmmsAttributeSet aAttrib;
      aAttrib.add(pwrs, "max");
      aAttrib.put(cur);
              
      shape="sphere";
      xmlNodePtr xmlCoefs = cur->xmlChildrenNode;
      while (xmlCoefs != NULL)
        {
          string cname((const char*)xmlCoefs->name);
          if (cname == "regions")
            {
              OhmmsAttributeSet cAttrib;
              cAttrib.add(shape, "shape");
              cAttrib.put(xmlCoefs);
              putContent(regions, xmlCoefs);
//               std::sort(regions.begin(),regions.end());
            }
          if (cname == "rho")  putContent(voffsets, xmlCoefs);
          xmlCoefs = xmlCoefs->next;
        }
        
      nreg=regions.size();
      app_log()<<" Shape: "<<shape<<endl;
      app_log()<<" Max Power: "<<pwrs<<endl;
      app_log()<<" regions:";
      for(int i(0);i<nreg;i++) app_log()<<" "<<regions[i];
      app_log()<<endl;
      
      
      if(voffsets.size()!=nreg*(grps+1))
        voffsets.resize(nreg*(grps+1),0);
      app_log()<<" voffsets:";
      for(int i(0);i<nreg*(grps+1);i++) app_log()<<" "<<voffsets[i];
      app_log()<<endl;
      
      Values.resize( (grps+1)*nreg*pwrs,0 );

      return true;
    }

    bool get(std::ostream& os) const {
      return true;
    }

    QMCHamiltonianBase* makeClone(ParticleSet& qp, TrialWaveFunction& psi)
    {
      NumberFluctuations* nf = new NumberFluctuations(qp);
      nf->nreg=nreg; nf->pwrs=pwrs; nf->grps=grps;
      nf->shape=shape; nf->myIndex=myIndex;
      nf->regions=regions; nf->voffsets=voffsets;
      nf->Values.resize( (grps+1)*nreg*pwrs,0 );
      
      return nf;
    }
    
  };
}
#endif

/***************************************************************************
 * $RCSfile$   $Author: jeongnim.kim $
 * $Revision: 4238 $   $Date: 2009-09-29 12:43:33 -0500 (Tue, 29 Sep 2009) $
 * $Id: CoulombPotential.h 4238 2009-09-29 17:43:33Z jeongnim.kim $ 
 ***************************************************************************/

