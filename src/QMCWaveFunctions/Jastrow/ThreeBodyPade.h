//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////
    
    
// BASED ON ONEBODYJASTROWFUNCTION.H; MODIFIED BY JOHN GERGELY

#ifndef QMCPLUSPLUS_GENERIC_THREEBODYPADE_H
#define QMCPLUSPLUS_GENERIC_THREEBODYPADE_H
#include "Configuration.h"
#include "QMCWaveFunctions/OrbitalBase.h"
#include "Particle/DistanceTableData.h"
#include "Particle/DistanceTable.h"
#include <fstream>
#include <iostream>

namespace qmcplusplus
{

//template<class FT>
class ThreeBodyPade: public OrbitalBase
{

  const DistanceTableData* d_table;
  const DistanceTableData* ee_table;

  struct Coeff
  {
    int m, n, o;
    ValueType c,d;
    Coeff(int m0, int n0, int o0, ValueType c0, ValueType d0): m(m0),n(n0),o(o0),c(c0),d(d0) {}
    //Coeff(){};
  };

  typedef int FT;
  ValueType curVal, curLap;
  GradType curGrad;
  ValueVectorType U,d2U;
  GradVectorType dU;
  ValueType *FirstAddressOfdU, *LastAddressOfdU;
  std::vector<FT*> Fs;
//    std::vector<FT*> Funique;
  std::vector<Coeff> C;
  // int N;  // number of terms in summation (user-defined)
public:

  typedef FT FuncType;

  inline void LoadC(int k, int m, int n, int o, double c)
  {
    C[k].m = m;
    C[k].n = n;
    C[k].o = o;
    C[k].c = c;
  }

  inline void InitC()
  {
    LoadC(0,0,0,1,0.5);
    LoadC(1,0,0,2,-.50464);
    LoadC(2,0,0,3,-.69909);
    LoadC(3,0,0,4,1.31627);
    LoadC(4,2,0,0,0.02031);
    LoadC(5,3,0,0,0.09953);
    LoadC(6,4,0,0,0.15156);
    LoadC(7,2,2,0,-3.01432);
    LoadC(8,2,0,2,1.44578);
    LoadC(9,2,2,2,7.26368);
    LoadC(10,4,0,2,-2.22981);
    LoadC(11,2,0,4,-2.52090);
    LoadC(12,4,2,2,-4.29786);
    LoadC(13,6,0,2,0.10429);
    LoadC(14,4,0,4,4.33268);
    LoadC(15,2,2,4,1.52266);
    LoadC(16,2,0,6,-1.15585);
    std::cerr << "Hard-wired coefficients: " << std::endl;
    std::cerr << "index m n o c" << std::endl;
    for(int k=0; k<C.size(); k++)
      std::cerr << k << "     " << C[k].m << " " << C[k].n << " " << C[k].o << " " << C[k].c << std::endl;
  }


  ///constructor
  ThreeBodyPade(ParticleSet& ions, ParticleSet& els) :  FirstAddressOfdU(NULL), LastAddressOfdU(NULL)
  {
    U.resize(els.getTotalNum());
    ee_table = DistanceTable::add(els);
    d_table = DistanceTable::add(ions,els);
    C.reserve(32); // need to reserve the memory...
    //InitC();
  }

  ~ThreeBodyPade() { }

  //evaluate the distance table with P
  void resetTargetParticleSet(ParticleSet& P)
  {
    d_table = DistanceTable::add(d_table->origin(),P);
  }

  void resetParameters(OptimizableSetType& optVariables)
  {
    ///Do nothing
  }

  inline ValueType Delta(int m, int n)
  {
    double value = 1;
    if(m == n)
      value = 0.5;
    return value;
  }

  ValueType evaluateLog(ParticleSet& P, ParticleSet::ParticleGradient_t& G, ParticleSet::ParticleLaplacian_t& L)
  {
    LogValue=0.0;
    U=0.0;
    int nnI,nnJ,foundI,foundJ;
    double rij, ri, rj, Rij, Ri, Rj, riInv, rjInv, rijInv; //to store distances and corresponding Pade terms
    for(int i=0; i<ee_table->size(SourceIndex); i++)
    {
      for(int ee=ee_table->M[i]; ee<ee_table->M[i+1]; ee++)
      {
        int j = ee_table->J[ee];
        nnI = i;
        nnJ = j;
        int powMax = 6;
        std::vector<double> RiPow(powMax);
        std::vector<double> RjPow(powMax);
        std::vector<double> RijPow(powMax);
        rij = ee_table->r(ee);
        ri = d_table->r(nnI);
        rj = d_table->r(nnJ);
        riInv = 1.0/(1+ri);
        rjInv = 1.0/(1+rj);
        rijInv = 1.0/(1+rij);
        Ri = ri*riInv;
        Rj = rj*rjInv;
        Rij = rij*rijInv;
//cerr << "ri = " << d_table->dr(nnI) << std::endl;
//cerr << "rj = " << d_table->dr(nnJ) << std::endl;
//cerr << "rij = " << ee_table->dr(ee) << std::endl; // store powers of the Pade terms rather than recompute
        for(int p = 0; p<=powMax; p++)
        {
          RiPow[p] = pow(Ri,p);
          RjPow[p] = pow(Rj,p);
          RijPow[p] = pow(Rij,p);
        }
        ValueType uij = 0;
        double lapI = 0;  // accumulate gradient and laplacian contributions
        double lapJ = 0;
        double grI = 0;
        double grJ = 0;
        double grIJ = 0;
        PosType drI = d_table->dr(nnI)*d_table->rinv(nnI); //unit vector
        PosType drJ = d_table->dr(nnJ)*d_table->rinv(nnJ);
        PosType drIJ = ee_table->dr(ee)*ee_table->rinv(ee); //unit vector
        int m,n,o;
        double c, del;
        // Loop over the N sets of coefficients for the correlation function
//
        for(int k = 0; k<C.size(); k++)
        {
          m = C[k].m;
          n = C[k].n;
          o = C[k].o;
          c = C[k].c;
          del = C[k].d;
//m = 0; n = 0; o=1; c=-1.0;
          //del = Delta(m,n); // =0.5, which compensates correctly for the extra factors of 2
          uij += del*c*(RiPow[m]*RjPow[n] + RjPow[m]*RiPow[n])*RijPow[o];
//
          grI += del*c*(m*RiPow[m-1]*RjPow[n] + n*RjPow[m]*RiPow[n-1])*RijPow[o]*riInv*riInv;
          grJ += del*c*(m*RjPow[m-1]*RiPow[n] + n*RiPow[m]*RjPow[n-1])*RijPow[o]*rjInv*rjInv;
          grIJ -= del*c*(RiPow[m]*RjPow[n] + RjPow[m]*RiPow[n])*(o*RijPow[o-1]*rijInv*rijInv);
          lapI += del*c*((m*(m-1)*RiPow[m-2]*RjPow[n] + n*(n-1)*RjPow[m]*RiPow[n-2])*RijPow[o]*riInv*riInv*riInv*riInv
                         + 2*(m*RiPow[m-1]*RjPow[n] + n*RjPow[m]*RiPow[n-1])*RijPow[o]*riInv*riInv*riInv*d_table->rinv(nnI)
                         - 2*(m*RiPow[m-1]*RjPow[n] + n*RjPow[m]*RiPow[n-1])*(o*RijPow[o-1])*(dot(drI,drIJ)*rijInv*rijInv*riInv*riInv)
                         + (RiPow[m]*RjPow[n] + RjPow[m]*RiPow[n])*(o*(o-1)*RijPow[o-2]*rijInv*rijInv*rijInv*rijInv)
                         + 2*o*(RiPow[m]*RjPow[n] + RjPow[m]*RiPow[n])*RijPow[o-1]*rijInv*rijInv*rijInv*ee_table->rinv(ee));
          lapJ += del*c*((m*(m-1)*RjPow[m-2]*RiPow[n] + n*(n-1)*RiPow[m]*RjPow[n-2])*RijPow[o]*rjInv*rjInv*rjInv*rjInv
                         + 2*(m*RjPow[m-1]*RiPow[n] + n*RiPow[m]*RjPow[n-1])*RijPow[o]*rjInv*rjInv*rjInv*d_table->rinv(nnJ)
                         + 2*(m*RjPow[m-1]*RiPow[n] + n*RiPow[m]*RjPow[n-1])*(o*RijPow[o-1])*(dot(drJ,drIJ)*rijInv*rijInv*rjInv*rjInv)
                         + (RiPow[m]*RjPow[n] + RjPow[m]*RiPow[n])*(o*(o-1)*RijPow[o-2]*rijInv*rijInv*rijInv*rijInv)
                         + 2*o*(RiPow[m]*RjPow[n] + RjPow[m]*RiPow[n])*RijPow[o-1]*rijInv*rijInv*rijInv*ee_table->rinv(ee));
//cerr << "lapI is now " << lapI << " and lapJ " << lapJ << std::endl;
//cerr << "  and dot test: " << d_table->dr(nnI) << "*" << ee_table->dr(ee) << " = " << dot(d_table->dr(nnI),ee_table->dr(ee)) << std::endl;
        }
        /*
        //////////////////////////////////////
        // just a chunk from TwoBody for testing
                  double dudr, d2udr2;
                  dudr = rijInv*rijInv;  // something like this...
                  d2udr2 = -2*dudr*rijInv;  // and this...
                  double testuij = -Rij;
                  dudr *= ee_table->rinv(ee);
                  PosType testdrI = ee_table->dr(ee);
                  PosType testdrJ = -1*drI;
                  double testgrI = dudr;
                  double testgrJ = grI;
                  double testlapI = -2.0*rijInv*rijInv*rijInv*ee_table->rinv(ee);
                  //lapI = d2udr2+2.0*dudr;
                  double testlapJ = lapI;
        ///////////////////////////////////
        */
        LogValue += uij;
        U[i] += uij;
        U[j] += uij;
        G[i] += drI*grI + drIJ*grIJ;
        G[j] += drJ*grJ - drIJ*grIJ;
        L[i] += lapI;
        L[j] += lapJ;
      }
    }
    return LogValue;
  }

  ValueType evaluate(ParticleSet& P,
                     ParticleSet::ParticleGradient_t& G,
                     ParticleSet::ParticleLaplacian_t& L)
  {
    return exp(evaluateLog(P,G,L));
  }

  /** evaluate the ratio \f$exp(U(iat)-U_0(iat))\f$
   * @param P active particle set
   * @param iat particle that has been moved.
   */
  inline ValueType ratio(ParticleSet& P, int iat)
  {
    return 0;
  }

  inline void restore(int iat) {}

// copied from OneBody... justn need it for now I think
  void addFunc(int source_type, FT* afunc)
  {
    const ParticleSet& ions=d_table->origin();
    if(Fs.empty())
    {
      Fs.resize(ions.getTotalNum(),0);
    }
    for(int i=0; i<Fs.size(); i++)
    {
      if(ions.GroupID[i] == source_type)
        Fs[i]=afunc;
    }
    //Funique.push_back(afunc);
  }


  /** equivalent to evalaute with additional data management */
  void registerData(ParticleSet& P, WFBufferType& buf)
  {
    //U.resize(d_table->size(VisitorIndex));
    d2U.resize(d_table->size(VisitorIndex));
    dU.resize(d_table->size(VisitorIndex));
    LogValue = 0.0;
    U=0.0;
    dU=0.0;
    d2U=0.0;
    ValueType uij, dudr, d2udr2;
    for(int i=0; i<d_table->size(SourceIndex); i++)
    {
      FT* func=Fs[i];
      if(func == 0)
        continue;
      for(int nn=d_table->M[i]; nn<d_table->M[i+1]; nn++)
      {
        int j = d_table->J[nn];
        //uij = F[d_table->PairID[nn]]->evaluate(d_table->r(nn), dudr, d2udr2);
        //uij = func->evaluate(d_table->r(nn), dudr, d2udr2);
        LogValue-=uij;
        U[j]+=uij;
        dudr *= d_table->rinv(nn);
        dU[j] -= dudr*d_table->dr(nn);
        d2U[j] -= d2udr2+2.0*dudr;
        //add gradient and laplacian contribution
        P.G[j] -= dudr*d_table->dr(nn);
        P.L[j] -= d2udr2+2.0*dudr;
      }
    }
    FirstAddressOfdU = &(dU[0][0]);
    LastAddressOfdU = FirstAddressOfdU + dU.size()*DIM;
    //add U, d2U and dU. Keep the order!!!
    buf.add(U.begin(), U.end());
    buf.add(d2U.begin(), d2U.end());
    buf.add(FirstAddressOfdU,LastAddressOfdU);
  }

  ValueType updateBuffer(ParticleSet& P, WFBufferType& buf)
  {
    LogValue = 0.0;
    U=0.0;
    dU=0.0;
    d2U=0.0;
    ValueType uij, dudr, d2udr2;
    for(int i=0; i<d_table->size(SourceIndex); i++)
    {
      FT* func=Fs[i];
      if(func == 0)
        continue;
      for(int nn=d_table->M[i]; nn<d_table->M[i+1]; nn++)
      {
        int j = d_table->J[nn];
        //uij = F[d_table->PairID[nn]]->evaluate(d_table->r(nn), dudr, d2udr2);
        //uij = func->evaluate(d_table->r(nn), dudr, d2udr2);
        LogValue-=uij;
        U[j]+=uij;
        dudr *= d_table->rinv(nn);
        dU[j] -= dudr*d_table->dr(nn);
        d2U[j] -= d2udr2+2.0*dudr;
        //add gradient and laplacian contribution
        P.G[j] -= dudr*d_table->dr(nn);
        P.L[j] -= d2udr2+2.0*dudr;
      }
    }
    FirstAddressOfdU = &(dU[0][0]);
    LastAddressOfdU = FirstAddressOfdU + dU.size()*DIM;
    buf.put(U.begin(), U.end());
    buf.put(d2U.begin(), d2U.end());
    buf.put(FirstAddressOfdU,LastAddressOfdU);
    return LogValue;
  }

  /** copy the current data from a buffer
   *@param P the ParticleSet to operate on
   *@param buf PooledData which stores the data for each walker
   *
   *copyFromBuffer uses the data stored by registerData or evaluate(P,buf)
   */
  void copyFromBuffer(ParticleSet& P, WFBufferType& buf)
  {
    buf.get(U.begin(), U.end());
    buf.get(d2U.begin(), d2U.end());
    buf.get(FirstAddressOfdU,LastAddressOfdU);
  }

  /** return the current value and copy the current data to a buffer
   *@param P the ParticleSet to operate on
   *@param buf PooledData which stores the data for each walker
   */
  inline ValueType evaluate(ParticleSet& P, PooledData<RealType>& buf)
  {
    ValueType sumu = 0.0;
    for(int i=0; i<U.size(); i++)
      sumu+=U[i];
    buf.put(U.begin(), U.end());
    buf.put(d2U.begin(), d2U.end());
    buf.put(FirstAddressOfdU,LastAddressOfdU);
    return exp(-sumu);
  }

  void acceptMove(ParticleSet& P, int iat)
  {
    std::cerr << "ThreeBodyPade::acceptMove does nothing right now" << std::endl;
  }

  void put(xmlNodePtr cur, VarRegistry<RealType>& vlist)
  {
    int m,n,o;
    ValueType c,d;
    const xmlChar* tpr = xmlGetProp(cur, (const xmlChar *)"id");
    const char* id=(const char*)tpr;
    tpr=xmlGetProp(cur, (const xmlChar *)"m");
    if(tpr != NULL)
    {
      m= atoi((const char*)tpr);
    }
    else
    {
      std::cerr << "PARSE ERROR: m not found." << std::endl;
    }
    tpr=xmlGetProp(cur, (const xmlChar *)"n");
    if(tpr != NULL)
    {
      n= atoi((const char*)tpr);
    }
    else
    {
      std::cerr << "PARSE ERROR: n not found." << std::endl;
    }
    tpr=xmlGetProp(cur, (const xmlChar *)"o");
    if(tpr != NULL)
    {
      o= atoi((const char*)tpr);
    }
    else
    {
      std::cerr << "PARSE ERROR: o not found." << std::endl;
    }
    tpr=xmlGetProp(cur, (const xmlChar *)"c");
    if(tpr != NULL)
    {
      c= atof((const char*)tpr);
    }
    else
    {
      std::cerr << "PARSE ERROR: c not found." << std::endl;
    }
    d = Delta(m,n);
    int now=C.size();
    Coeff newCoeff(m,n,o,c,d);
    C.push_back(newCoeff);
    //sprintf(ida.c_str(),id);
    std::cerr << "Adding to VarList: " << (const char*)id << ", " << c << " located at " << &C[C.size()].c << " and index " << C.size() << std::endl;
    vlist[(const char*)id]=C[now].c;
    //vlist.add((const char*)id,&(C[now].c),1);
    XMLReport("Pade (A*r+C*r*r)/(1+Br) = (" << A << "," << B << "," << C << ")")
    std::cerr << "Added coefficient set " << newCoeff.m << ", " << newCoeff.n << ", " << newCoeff.o << ", " << newCoeff.c << ", " << newCoeff.d  << std::endl;
  }
};
}
#endif
