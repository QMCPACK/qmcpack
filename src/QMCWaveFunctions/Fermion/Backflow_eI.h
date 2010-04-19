//////////////////////////////////////////////////////////////////
// (c) Copyright 2003-  by Jeongnim Kim
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
#ifndef QMCPLUSPLUS_BACKFLOW_ELEC-ION_H
#define QMCPLUSPLUS_BACKFLOW_ELEC-ION_H
#include "QMCWaveFunctions/OrbitalSetTraits.h"
#include "QMCWaveFunctions/Fermion/BackflowFunctionBase.h"
#include "Particle/DistanceTable.h"
#include <cmath>
#include <vector>

namespace qmcplusplus
{
  template<class FT>
  class Backflow_eI: public BackflowFunctionBase 
  {

    public:

    bool uniqueFunctions;
 
    vector<FT*> RadFun;

    Backflow_eI(ParticleSet& ions, ParticleSet& els): BackflowFunctionBase(ions,els),uniqueFunctions(false) 
    {
      myTable = DistanceTable::add(ions,els); 
    }

    Backflow_eI(ParticleSet& ions, ParticleSet& els, FT* RF): BackflowFunctionBase(ions,els),uniqueFunctions(false)
    {
      myTable = DistanceTable::add(ions,els);
      // same radial function for all centers
      for(int i=0; i<NumCenters; i++) RadFun.push_back(RF);
      uniqueFunctions=false;
    }

    ~Backflow_eI() {}; 
 
    void resetParameters(const opt_variables_type& active)
    {
       RadFun[0]->resetParameters(active);
       if(uniqueFunctions)
        for(int i=1; i<RadFun.size(); i++) RadFun[i]->resetParameters(active);
    }

    void checkInVariables(opt_variables_type& active)
    {
      RadFun[0]->checkInVariables(active);
       if(uniqueFunctions)
        for(int i=1; i<RadFun.size(); i++) RadFun[i]->checkInVariables(active);
    }

    void checkOutVariables(const opt_variables_type& active)
    {
      RadFun[0]->checkOutVariables(active);
       if(uniqueFunctions)
        for(int i=1; i<RadFun.size(); i++) RadFun[i]->checkOutVariables(active);
    }


    BackflowFunctionBase* makeClone()
    {
       Backflow_eI* clone = new Backflow_eI(*this);
       return clone; 
    }

    
    /** calculate quasi-particle coordinates only
     */
    inline void 
    evaluate(const ParticleSet& P, ParticleSet& QP)
    {
      ValueType du,d2u; 
      for(int i=0; i<myTable->size(SourceIndex); i++) {
        for(int nn=myTable->M[i]; nn<myTable->M[i+1]; nn++) {
          int j = myTable->J[nn];
          ValueType uij = RadFun[i]->evaluate(myTable->r(nn),du,d2u);
          PosType u = uij*myTable->dr(nn);
          QP.R[j] += u;  // dr(ij) = r_j-r_i 
        }
      }
    }


    inline void
    evaluate(const ParticleSet& P, ParticleSet& QP, GradVector_t& Bmat, HessMatrix_t& Amat)
    {
      ValueType du,d2u,temp;
      for(int i=0; i<myTable->size(SourceIndex); i++) {
        for(int nn=myTable->M[i]; nn<myTable->M[i+1]; nn++) {
          int j = myTable->J[nn];
          ValueType uij = RadFun[i]->evaluate(myTable->r(nn),du,d2u);
          PosType u = uij*myTable->dr(nn);
          du *= myTable->rinv(nn);
          QP.R[j] += u;

          RealType *ptj = Amat(j,j).begin();
          temp=du*myTable->dr(nn)[0]*myTable->dr(nn)[0]+uij;  //(0,0)
          *(ptj) += temp;
          temp=du*myTable->dr(nn)[1]*myTable->dr(nn)[1]+uij;  //(1,1)
          *(ptj+4) += temp;
          temp=du*myTable->dr(nn)[2]*myTable->dr(nn)[2]+uij;  //(2,2)
          *(ptj+8) += temp;
          temp=du*myTable->dr(nn)[1]*myTable->dr(nn)[0];  //(1,0)+(0,1)
          *(ptj+1) += temp;
          *(ptj+3) += temp;
          temp=du*myTable->dr(nn)[2]*myTable->dr(nn)[0];  //(2,0)+(2,0)
          *(ptj+2) += temp;
          *(ptj+6) += temp;
          temp=du*myTable->dr(nn)[1]*myTable->dr(nn)[2];  //(1,2)+(2,1)
          *(ptj+5) += temp;
          *(ptj+7) += temp;

          u = (d2u+4.0*du)*myTable->dr(nn);
          Bmat(j) += u;
        }
      }
    }

    
    /** calculate quasi-particle coordinates, Bmat and Amat 
     */
    inline void
    evaluate(const ParticleSet& P, ParticleSet& QP, GradMatrix_t& Bmat_full, HessMatrix_t& Amat)
    {
      RealType du,d2u,temp;
      for(int i=0; i<myTable->size(SourceIndex); i++) {
        for(int nn=myTable->M[i]; nn<myTable->M[i+1]; nn++) {
          int j = myTable->J[nn];
          ValueType uij = RadFun[i]->evaluate(myTable->r(nn),du,d2u);
          du *= myTable->rinv(nn);
          PosType u = uij*myTable->dr(nn); 
          QP.R[j] += u; 
 
          RealType *ptj = Amat(j,j).begin(); 
          temp=du*myTable->dr(nn)[0]*myTable->dr(nn)[0]+uij;  //(0,0)
          *(ptj) += temp;
          temp=du*myTable->dr(nn)[1]*myTable->dr(nn)[1]+uij;  //(1,1)
          *(ptj+4) += temp;
          temp=du*myTable->dr(nn)[2]*myTable->dr(nn)[2]+uij;  //(2,2)
          *(ptj+8) += temp;
          temp=du*myTable->dr(nn)[1]*myTable->dr(nn)[0];  //(1,0)+(0,1)
          *(ptj+1) += temp;
          *(ptj+3) += temp;
          temp=du*myTable->dr(nn)[2]*myTable->dr(nn)[0];  //(2,0)+(2,0)
          *(ptj+2) += temp;
          *(ptj+6) += temp;
          temp=du*myTable->dr(nn)[1]*myTable->dr(nn)[2];  //(1,2)+(2,1)
          *(ptj+5) += temp;
          *(ptj+7) += temp;

// this will create problems with QMC_COMPLEX, because Bmat is ValueType and dr is RealType
          u = (d2u+4.0*du)*myTable->dr(nn); 
          Bmat_full(j,j) += u;  
        }
      }
    }
  };

}

#endif
