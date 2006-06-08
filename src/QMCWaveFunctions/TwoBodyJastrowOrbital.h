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
#ifndef QMCPLUSPLUS_COMMON_TWOBODYJASTROW_H
#define QMCPLUSPLUS_COMMON_TWOBODYJASTROW_H
#include "Configuration.h"
#include  <map>
#include "QMCWaveFunctions/OrbitalBase.h"
#include <iostream>
#include <fstream>

namespace qmcplusplus {

  /** @ingroup OrbitalComponent
   *  @brief Specialization for two-body Jastrow function using multiple functors
   *
   *Each pair-type can have distinct function \f$u(r_{ij})\f$.
   *For electrons, distinct pair correlation functions are used 
   *for spins up-up/down-down and up-down/down-up.
   */ 
  template<class FT>
  class TwoBodyJastrowOrbital: public OrbitalBase {

    const DistanceTableData* d_table;

    int N,NN;
    ValueType DiffVal, DiffValSum;
    ValueVectorType U,d2U,curLap,curVal;
    GradVectorType dU,curGrad;
    ValueType *FirstAddressOfdU, *LastAddressOfdU;
    Matrix<int> PairID;

    std::map<std::string,FT*> J2Unique;

  public:

    typedef FT FuncType;
    //QIO out;
    ///container for the Jastrow functions 
    vector<FT*> F;
    ofstream out;
    ///constructor
    TwoBodyJastrowOrbital(ParticleSet& p, DistanceTableData* dtable): d_table(dtable) { 
      N=p.getTotalNum();
      NN=N*N;
      U.resize(NN+1);
      PairID.resize(N,N);
      int nsp=p.groups();
      for(int i=0; i<N; i++)
	for(int j=0; j<N; j++) 
	  PairID(i,j) = p.GroupID[i]*nsp+p.GroupID[j];
      out.open("TwoBodJast.dat");
      //out.Open("TwoBodJast.dat");
    }

    ~TwoBodyJastrowOrbital(){ }

    void insert(const string& aname, FT* j) {
      J2Unique[aname]=j;
    }

    void insert(std::map<std::string,FT*>& j2unique) {
      J2Unique.insert(j2unique.begin(),j2unique.end());
    }

    ///reset the value of all the unique Two-Body Jastrow functions
    void reset() { 
      typename std::map<std::string,FT*>::iterator it(J2Unique.begin()),it_end(J2Unique.end());
      while(it != it_end) {
        (*it).second->reset(); ++it;
      }
    }

    //Add a functor
    void addFunc(FT* afunc) { F.push_back(afunc); }

    //evaluate the distance table with els
    void resetTargetParticleSet(ParticleSet& P) {
      d_table = DistanceTable::getTable(DistanceTable::add(P));
    }

    /** 
     *@param P input configuration containing N particles
     *@param G a vector containing N gradients
     *@param L a vector containing N laplacians
     *@param G returns the gradient \f$G[i]={\bf \nabla}_i J({\bf R})\f$
     *@param L returns the laplacian \f$L[i]=\nabla^2_i J({\bf R})\f$
     *@return \f$exp(-J({\bf R}))\f$
     *@brief While evaluating the value of the Jastrow for a set of
     *particles add the gradient and laplacian contribution of the
     *Jastrow to G(radient) and L(aplacian) for local energy calculations
     *such that \f[ G[i]+={\bf \nabla}_i J({\bf R}) \f] 
     *and \f[ L[i]+=\nabla^2_i J({\bf R}). \f]
     *@note The DistanceTableData contains only distinct pairs of the 
     *particles belonging to one set, e.g., SymmetricDTD.
     */
    ValueType evaluateLog(ParticleSet& P,
		          ParticleSet::ParticleGradient_t& G, 
		          ParticleSet::ParticleLaplacian_t& L) {
      LogValue=0.0;
      ValueType dudr, d2udr2;
      PosType gr;
      for(int i=0; i<d_table->size(SourceIndex); i++) {
	for(int nn=d_table->M[i]; nn<d_table->M[i+1]; nn++) {
	  int j = d_table->J[nn];
	  //LogValue -= F[d_table->PairID[nn]]->evaluate(d_table->r(nn), dudr, d2udr2);
	  ValueType uij = F[d_table->PairID[nn]]->evaluate(d_table->r(nn), dudr, d2udr2);
          LogValue -= uij;
	  U[i*N+j]=uij; U[j*N+i]=uij; //save for the ratio
	  //multiply 1/r
	  dudr *= d_table->rinv(nn);
	  gr = dudr*d_table->dr(nn);
	  //(d^2 u \over dr^2) + (2.0\over r)(du\over\dr)
	  ValueType lap = d2udr2+2.0*dudr;

	  //multiply -1
	  G[i] += gr;
	  G[j] -= gr;
	  L[i] -= lap; 
	  L[j] -= lap; 
	}
      }
      //out.WriteLine(LogValue);
      //out << LogValue << endl;
      return LogValue;
    }

    ValueType evaluate(ParticleSet& P,
		       ParticleSet::ParticleGradient_t& G, 
		       ParticleSet::ParticleLaplacian_t& L) {
      return std::exp(evaluateLog(P,G,L));
    }

    ValueType ratio(ParticleSet& P, int iat) {
      ValueType d(0.0);
      const int* pairid(PairID[iat]);
      for(int jat=0, ij=iat*N; jat<N; jat++,ij++) {
        if(iat != jat) {
          d += U[ij]-F[pairid[jat]]->evaluate(d_table->Temp[jat].r1);
        }
      }
      return std::exp(d);
    }

    /** later merge the loop */
    ValueType ratio(ParticleSet& P, int iat,
		    ParticleSet::ParticleGradient_t& dG,
		    ParticleSet::ParticleLaplacian_t& dL)  {
      register ValueType dudr, d2udr2,u;
      register PosType gr;
      DiffVal = 0.0;      
      const int* pairid = PairID[iat];
      for(int jat=0, ij=iat*N; jat<N; jat++,ij++) {
	if(jat==iat) {
	  curVal[jat] = 0.0;curGrad[jat]=0.0; curLap[jat]=0.0;
	} else {
	  curVal[jat] = F[pairid[jat]]->evaluate(d_table->Temp[jat].r1, dudr, d2udr2);
	  dudr *= d_table->Temp[jat].rinv1;
	  curGrad[jat] = -dudr*d_table->Temp[jat].dr1;
	  curLap[jat] = -(d2udr2+2.0*dudr);
	  DiffVal += (U[ij]-curVal[jat]);
	}
      }
      GradType sumg,dg;
      ValueType suml=0.0,dl;
      for(int jat=0,ij=iat*N,ji=iat; jat<N; jat++,ij++,ji+=N) {
	sumg += (dg=curGrad[jat]-dU[ij]);
	suml += (dl=curLap[jat]-d2U[ij]);
        dG[jat] -= dg;
	dL[jat] += dl;
      }
      dG[iat] += sumg;
      dL[iat] += suml;     
      return std::exp(DiffVal);
    }

    /** later merge the loop */
    ValueType logRatio(ParticleSet& P, int iat,
		    ParticleSet::ParticleGradient_t& dG,
		    ParticleSet::ParticleLaplacian_t& dL)  {
      register ValueType dudr, d2udr2,u;
      register PosType gr;
      DiffVal = 0.0;      
      const int* pairid = PairID[iat];
      for(int jat=0, ij=iat*N; jat<N; jat++,ij++) {
	if(jat==iat) {
	  curVal[jat] = 0.0;curGrad[jat]=0.0; curLap[jat]=0.0;
	} else {
	  curVal[jat] = F[pairid[jat]]->evaluate(d_table->Temp[jat].r1, dudr, d2udr2);
	  dudr *= d_table->Temp[jat].rinv1;
	  curGrad[jat] = -dudr*d_table->Temp[jat].dr1;
	  curLap[jat] = -(d2udr2+2.0*dudr);
	  DiffVal += (U[ij]-curVal[jat]);
	}
      }
      GradType sumg,dg;
      ValueType suml=0.0,dl;
      for(int jat=0,ij=iat*N,ji=iat; jat<N; jat++,ij++,ji+=N) {
	sumg += (dg=curGrad[jat]-dU[ij]);
	suml += (dl=curLap[jat]-d2U[ij]);
        dG[jat] -= dg;
	dL[jat] += dl;
      }
      dG[iat] += sumg;
      dL[iat] += suml;     
      return DiffVal;
    }

    inline void restore(int iat) {}

    void acceptMove(ParticleSet& P, int iat) { 
      DiffValSum += DiffVal;
      for(int jat=0,ij=iat*N,ji=iat; jat<N; jat++,ij++,ji+=N) {
	dU[ij]=curGrad[jat]; 
	dU[ji]=-1.0*curGrad[jat];
	d2U[ij]=d2U[ji] = curLap[jat];
	U[ij] =  U[ji] = curVal[jat];
      }
    }


    inline void update(ParticleSet& P, 		
		       ParticleSet::ParticleGradient_t& dG, 
		       ParticleSet::ParticleLaplacian_t& dL,
		       int iat) {
      DiffValSum += DiffVal;
      GradType sumg,dg;
      ValueType suml=0.0,dl;
      for(int jat=0,ij=iat*N,ji=iat; jat<N; jat++,ij++,ji+=N) {
	sumg += (dg=curGrad[jat]-dU[ij]);
	suml += (dl=curLap[jat]-d2U[ij]);
	dU[ij]=curGrad[jat]; 
	dU[ji]=-1.0*curGrad[jat];
	d2U[ij]=d2U[ji] = curLap[jat];
	U[ij] =  U[ji] = curVal[jat];
        dG[jat] -= dg;
	dL[jat] += dl;
      }
      dG[iat] += sumg;
      dL[iat] += suml;     
    }

    inline ValueType registerData(ParticleSet& P, PooledData<RealType>& buf){
      N=d_table->size(VisitorIndex);
      NN=N*N;
      U.resize(NN+1);
      d2U.resize(NN);
      dU.resize(NN);
      curGrad.resize(N);
      curLap.resize(N);
      curVal.resize(N);
      ValueType dudr, d2udr2,u;

      LogValue=0.0;
      PosType gr;
      PairID.resize(d_table->size(SourceIndex),d_table->size(SourceIndex));
      int nsp=P.groups();
      for(int i=0; i<d_table->size(SourceIndex); i++)
	for(int j=0; j<d_table->size(SourceIndex); j++) 
	  PairID(i,j) = P.GroupID[i]*nsp+P.GroupID[j];

      for(int i=0; i<d_table->size(SourceIndex); i++) {
	for(int nn=d_table->M[i]; nn<d_table->M[i+1]; nn++) {
	  int j = d_table->J[nn];
	  //ValueType sumu = F.evaluate(d_table->r(nn));
	  u = F[d_table->PairID[nn]]->evaluate(d_table->r(nn), dudr, d2udr2);
	  LogValue -= u;
	  dudr *= d_table->rinv(nn);
	  gr = dudr*d_table->dr(nn);
	  //(d^2 u \over dr^2) + (2.0\over r)(du\over\dr)\f$
	  ValueType lap = d2udr2+2.0*dudr;
	  int ij = i*N+j, ji=j*N+i;
	  U[ij]=u; U[ji]=u;
	  dU[ij] = gr; dU[ji] = -1.0*gr;
	  d2U[ij] = -lap; d2U[ji] = -lap;

	  //add gradient and laplacian contribution
	  P.G[i] += gr;
	  P.G[j] -= gr;
	  P.L[i] -= lap; 
	  P.L[j] -= lap; 
	}
      }

      U[NN]= LogValue;
      buf.add(U.begin(), U.end());
      buf.add(d2U.begin(), d2U.end());
      FirstAddressOfdU = &(dU[0][0]);
      LastAddressOfdU = FirstAddressOfdU + dU.size()*DIM;
      buf.add(FirstAddressOfdU,LastAddressOfdU);

      return LogValue;
    }

    inline ValueType updateBuffer(ParticleSet& P, PooledData<RealType>& buf){
      ValueType dudr, d2udr2,u;
      LogValue=0.0;
      PosType gr;
      PairID.resize(d_table->size(SourceIndex),d_table->size(SourceIndex));
      int nsp=P.groups();
      for(int i=0; i<d_table->size(SourceIndex); i++)
	for(int j=0; j<d_table->size(SourceIndex); j++) 
	  PairID(i,j) = P.GroupID[i]*nsp+P.GroupID[j];

      for(int i=0; i<d_table->size(SourceIndex); i++) {
	for(int nn=d_table->M[i]; nn<d_table->M[i+1]; nn++) {
	  int j = d_table->J[nn];
	  //ValueType sumu = F.evaluate(d_table->r(nn));
	  u = F[d_table->PairID[nn]]->evaluate(d_table->r(nn), dudr, d2udr2);
	  LogValue -= u;
	  dudr *= d_table->rinv(nn);
	  gr = dudr*d_table->dr(nn);
	  //(d^2 u \over dr^2) + (2.0\over r)(du\over\dr)\f$
	  ValueType lap = d2udr2+2.0*dudr;
	  int ij = i*N+j, ji=j*N+i;
	  U[ij]=u; U[ji]=u;
	  dU[ij] = gr; dU[ji] = -1.0*gr;
	  d2U[ij] = -lap; d2U[ji] = -lap;

	  //add gradient and laplacian contribution
	  P.G[i] += gr;
	  P.G[j] -= gr;
	  P.L[i] -= lap; 
	  P.L[j] -= lap; 
	}
      }

      U[NN]= LogValue;
      buf.put(U.begin(), U.end());
      buf.put(d2U.begin(), d2U.end());
      buf.put(FirstAddressOfdU,LastAddressOfdU);
      return LogValue;
    }
    
    inline void copyFromBuffer(ParticleSet& P, PooledData<RealType>& buf) {
      buf.get(U.begin(), U.end());
      buf.get(d2U.begin(), d2U.end());
      buf.get(FirstAddressOfdU,LastAddressOfdU);
      DiffValSum=0.0;
    }

    inline ValueType evaluate(ParticleSet& P, PooledData<RealType>& buf) {
      ValueType x = (U[NN] += DiffValSum);
      buf.put(U.begin(), U.end());
      buf.put(d2U.begin(), d2U.end());
      buf.put(FirstAddressOfdU,LastAddressOfdU);
      return std::exp(x); 
    }

  };
}
#endif
/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$ 
 ***************************************************************************/

