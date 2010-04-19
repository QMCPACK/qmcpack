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
#ifndef QMCPLUSPLUS_BACKFLOW_E-E_H
#define QMCPLUSPLUS_BACKFLOW_E-E_H
#include "QMCWaveFunctions/OrbitalSetTraits.h"
#include "QMCWaveFunctions/Fermion/BackflowFunctionBase.h"
#include "Particle/DistanceTable.h"
#include <cmath>

namespace qmcplusplus
{
  template<class FT>
  class Backflow_ee: public BackflowFunctionBase 
  {

    public:

    FT *RadFun;

    Backflow_ee(ParticleSet& ions, ParticleSet& els): BackflowFunctionBase(ions,els),RadFun(0) 
    {
      myTable = DistanceTable::add(els,els);
    }

    Backflow_ee(ParticleSet& ions, ParticleSet& els, FT* RF): BackflowFunctionBase(ions,els),RadFun(RF) 
    {
      myTable = DistanceTable::add(els,els);
    }

    ~Backflow_ee() {}; 
 
    void resetParameters(const opt_variables_type& active)
    {
       RadFun->resetParameters(active);
    }

    void checkInVariables(opt_variables_type& active)
    {
      RadFun->checkInVariables(active);
    }

    void checkOutVariables(const opt_variables_type& active)
    {
      RadFun->checkOutVariables(active);
    }


    BackflowFunctionBase* makeClone()
    {
       Backflow_ee* clone = new Backflow_ee(*this);
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
          ValueType uij = RadFun->evaluate(myTable->r(nn),du,d2u);
          PosType u = uij*myTable->dr(nn);
          QP.R[i] -= u;  // dr(ij) = r_j-r_i 
          QP.R[j] += u;  
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
          ValueType uij = RadFun->evaluate(myTable->r(nn),du,d2u);
          PosType u = uij*myTable->dr(nn);
          du *= myTable->rinv(nn);
          QP.R[i] -= u;
          QP.R[j] += u;

// assume 3D for now, add compiler directives later 
// access of Tensor is column based (f77 style)
// can I do an outer product of dr to get a Tensor?? Is it implemented??
          RealType *pt = Amat(i,j).begin();
          RealType *pti = Amat(i,i).begin();
          RealType *ptj = Amat(j,j).begin();
          temp=du*myTable->dr(nn)[0]*myTable->dr(nn)[0]+uij;  //(0,0)
          *(pt) = -temp;
          *(pti) += temp;
          *(ptj) += temp;
          temp=du*myTable->dr(nn)[1]*myTable->dr(nn)[1]+uij;  //(1,1)
          *(pt+4) = -temp;
          *(pti+4) += temp;
          *(ptj+4) += temp;
          temp=du*myTable->dr(nn)[2]*myTable->dr(nn)[2]+uij;  //(2,2)
          *(pt+8) = -temp;
          *(pti+8) += temp;
          *(ptj+8) += temp;
          temp=du*myTable->dr(nn)[1]*myTable->dr(nn)[0];  //(1,0)+(0,1)
          *(pt+1) = -temp;
          *(pt+3) = -temp;
          *(pti+1) += temp;
          *(pti+3) += temp;
          *(ptj+1) += temp;
          *(ptj+3) += temp;
          temp=du*myTable->dr(nn)[2]*myTable->dr(nn)[0];  //(2,0)+(2,0)
          *(pt+2) = -temp;
          *(pt+6) = -temp;
          *(pti+2) += temp;
          *(pti+6) += temp;
          *(ptj+2) += temp;
          *(ptj+6) += temp;
          temp=du*myTable->dr(nn)[1]*myTable->dr(nn)[2];  //(1,2)+(2,1)
          *(pt+5) = -temp;
          *(pt+7) = -temp;
          *(pti+5) += temp;
          *(pti+7) += temp;
          *(ptj+5) += temp;
          *(ptj+7) += temp;
          Amat(j,i)=Amat(i,j);  // symmetrize wrt i and j

// this will create problems with QMC_COMPLEX, because Bmat is ValueType and dr is RealType
          u = 2.0*(d2u+4.0*du)*myTable->dr(nn);
          Bmat(i) -= u;
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
          ValueType uij = RadFun->evaluate(myTable->r(nn),du,d2u);
          du *= myTable->rinv(nn);
          PosType u = uij*myTable->dr(nn); 
          QP.R[i] -= u;  
          QP.R[j] += u; 
 
// assume 3D for now, add compiler directives later 
          RealType *pt = Amat(i,j).begin(); 
          RealType *pti = Amat(i,i).begin(); 
          RealType *ptj = Amat(j,j).begin(); 
          temp=du*myTable->dr(nn)[0]*myTable->dr(nn)[0]+uij;  //(0,0)
          *(pt) = -temp;    
          *(pti) += temp;
          *(ptj) += temp;
          temp=du*myTable->dr(nn)[1]*myTable->dr(nn)[1]+uij;  //(1,1)
          *(pt+4) = -temp;    
          *(pti+4) += temp;
          *(ptj+4) += temp;
          temp=du*myTable->dr(nn)[2]*myTable->dr(nn)[2]+uij;  //(2,2)
          *(pt+8) = -temp;    
          *(pti+8) += temp;
          *(ptj+8) += temp;
          temp=du*myTable->dr(nn)[1]*myTable->dr(nn)[0];  //(1,0)+(0,1)
          *(pt+1) = -temp;   
          *(pt+3) = -temp;   
          *(pti+1) += temp;
          *(pti+3) += temp;
          *(ptj+1) += temp;
          *(ptj+3) += temp;
          temp=du*myTable->dr(nn)[2]*myTable->dr(nn)[0];  //(2,0)+(2,0)
          *(pt+2) = -temp;    
          *(pt+6) = -temp;
          *(pti+2) += temp;
          *(pti+6) += temp;
          *(ptj+2) += temp;
          *(ptj+6) += temp;
          temp=du*myTable->dr(nn)[1]*myTable->dr(nn)[2];  //(1,2)+(2,1)
          *(pt+5) = -temp;    
          *(pt+7) = -temp;
          *(pti+5) += temp;
          *(pti+7) += temp;
          *(ptj+5) += temp;
          *(ptj+7) += temp;
          Amat(j,i)=Amat(i,j);  // symmetrize wrt i and j

// this will create problems with QMC_COMPLEX, because Bmat is ValueType and dr is RealType
          // d2u + (ndim+1)*du
          u = (d2u+4.0*du)*myTable->dr(nn); 
          Bmat_full(i,i) -= u;  
          Bmat_full(j,j) += u;  
          Bmat_full(i,j) = u;  
          Bmat_full(j,i) = -1.0*u;
        }
      }
    }


  };

}

#endif
