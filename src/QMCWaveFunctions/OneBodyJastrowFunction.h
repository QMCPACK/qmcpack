//////////////////////////////////////////////////////////////////
// (c) Copyright 2003-  by Jeongnim Kim
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
#ifndef OHMMS_QMC_GENERIC_ONEBODYJASTROW_H
#define OHMMS_QMC_GENERIC_ONEBODYJASTROW_H
#include "Configuration.h"
#include "QMCWaveFunctions/OrbitalBase.h"

namespace ohmmsqmc {
  
  /** generic implementation of OneBodyJastrow<FT,SharedFunction>
   *
   *The One-Body Jastrow has the form  
   \f[ J_{e}({\bf R}) = \sum_{I=1}^{N_I} 
   \sum_{i=1}^{N_e} u(r_{iI}) \f]
   where \f[ r_{iI} = |{\bf r}_i - {\bf R}_I| \f]
   and the first summnation \f$ 1..N_I \f$ is over the centers
   while the second summnation \f$ 1..N_e \f$ is over electrons.
   *
   *The template parameter FT is a functional \f$ u(r_{iI}), \f$ e.g.
   *class PadeJastrow<T> or NoCuspJastrow<T>
   *Requriement of the template function is 
   *ValueType evaluate(ValueType r, ValueType& dudr, ValueType& d2udr2).
   *
   To calculate the Gradient use the identity
   \f[ {\bf \nabla}_i(r_{iI}) = \frac{{{\bf r_{iI}}}}{r_{iI}}. \f] 
   \f[ {\bf \nabla}_i(J_{e}({\bf R})= 
   \sum_{I=1}^{N_I} \frac{du}{dr_{iI}}{\bf \hat{r_{iI}}}.
   \f]
   To calculate The Laplacian use the identity
   \f[ \nabla^2_i(r_{iI})=\frac{2}{r_{iI}}, \f]
   and the vector product rule
   \f[ 
   \nabla^2 \cdot (f{\bf A}) =
   f({\bf \nabla} \cdot {\bf A})+{\bf A}\cdot({\bf \nabla} f)
   \f]
   \f[
   \nabla^2_i (J_{e}({\bf R})) = \sum_{I=1}^{N_I} 
   \left(\frac{du}{dr_{iI}}\right) {\bf \nabla}_i
   \cdot {\bf \hat{r_{iI}}} - {\bf \hat{r_{iI}}} \cdot 
   \left(\frac{d^2u}{dr_{iI}^2}\right){\bf \hat{r_{iI}}}
   \f]
   which can be simplified to
   \f[
   \nabla^2_i (J_{e}({\bf R})) = \sum_{I=1}^{N_I}
   \left(\frac{2}{r_{iI}}\frac{du}{dr_{iI}}
   + \frac{d^2u}{dr_{iI}^2}\right)
   \f]
   *
   *A generic OneBodyJastrow function uses a distance table.
   *The indices I(sources) and i(targets) are distinct. In general, the source 
   *particles are fixed, e.g., the nuclei, while the target particles are updated
   *by MC methods.
   *
   *@todo The second template parameter boolean will be removed, after
   *the efficieny of different vectorized schemes is evaluated.
   */
  template<class FT>
  class OneBodyJastrow: public OrbitalBase {

    const DistanceTableData* d_table;

    ValueType curVal, curLap;
    GradType curGrad;
    ValueVectorType U,d2U;
    GradVectorType dU;
    ValueType *FirstAddressOfdU, *LastAddressOfdU;

  public:

    typedef FT FuncType;

    vector<FT*> F;

    ///constructor
    OneBodyJastrow(DistanceTableData* dtable)
      : d_table(dtable), FirstAddressOfdU(NULL), LastAddressOfdU(NULL){ 

    }

    ~OneBodyJastrow(){
      DEBUGMSG("OneBodyJastrow::~OneBodyJastrow")
	//for(int i=0; i<F.size(); i++) delete F[i];
    }

    void reset() { 
      for(int i=0; i<F.size(); i++) F[i]->reset();
    }

    /** 
     *@param P input configuration containing N particles
     *@param G a vector containing N gradients
     *@param L a vector containing N laplacians
     *@return The wavefunction value  \f$exp(-J({\bf R}))\f$
     *
     *Upon exit, the gradient \f$G[i]={\bf \nabla}_i J({\bf R})\f$
     *and the laplacian \f$L[i]=\nabla^2_i J({\bf R})\f$ are accumulated.
     *While evaluating the value of the Jastrow for a set of
     *particles add the gradient and laplacian contribution of the
     *Jastrow to G(radient) and L(aplacian) for local energy calculations
     *such that \f[ G[i]+={\bf \nabla}_i J({\bf R}) \f] 
     *and \f[ L[i]+=\nabla^2_i J({\bf R}). \f]
     */
    ValueType evaluate(ParticleSet& P,
		       ParticleSet::ParticleGradient_t& G, 
		       ParticleSet::ParticleLaplacian_t& L) {
      ValueType sumu = 0.0;
      ValueType dudr, d2udr2;
      for(int i=0; i<d_table->size(SourceIndex); i++) {
	for(int nn=d_table->M[i]; nn<d_table->M[i+1]; nn++) {
	  int j = d_table->J[nn];
	  sumu += F[d_table->PairID[nn]]->evaluate(d_table->r(nn), dudr, d2udr2);
	  dudr *= d_table->rinv(nn);
	  G[j] -= dudr*d_table->dr(nn);
	  L[j] -= d2udr2+2.0*dudr;
	}
      }
      return exp(-sumu);
    }


    ValueType ratio(ParticleSet& P, int iat,
		    ParticleSet::ParticleGradient_t& dG,
		    ParticleSet::ParticleLaplacian_t& dL)  {
      ValueType v=ratio(P,iat);    
      dG[iat] += curGrad-dU[iat];
      dL[iat] += curLap-d2U[iat]; 
      return v;
    }

    inline void restore(int iat) {}

    void update(ParticleSet& P, int iat) {
      U[iat] = curVal;
      dU[iat]=curGrad;
      d2U[iat]=curLap;
    }
    

    ValueType ratio(ParticleSet& P, int iat) {
      int n=d_table->size(VisitorIndex);
      curVal=0.0;
      curLap=0.0;
      curGrad = 0.0;
      ValueType dudr, d2udr2;
      for(int i=0, nn=iat; i<d_table->size(SourceIndex); i++,nn+= n) {
	int ij=d_table->PairID[nn];
	curVal += F[ij]->evaluate(d_table->Temp[i].r1,dudr,d2udr2);
	dudr *= d_table->Temp[i].rinv1;
	curGrad -= dudr*d_table->Temp[i].dr1;
	curLap  -= d2udr2+2.0*dudr;
      }
      return exp(U[iat]-curVal);
    } 	  

    void update(ParticleSet& P, 
		ParticleSet::ParticleGradient_t& dG, 
		ParticleSet::ParticleLaplacian_t& dL,
		int iat) {
      dG[iat] += curGrad-dU[iat]; dU[iat]=curGrad;
      dL[iat] += curLap-d2U[iat]; d2U[iat]=curLap;
      U[iat] = curVal;
    }


    /** equivalent to evalaute with additional data management */
    void registerData(ParticleSet& P, PooledData<RealType>& buf){

      U.resize(d_table->size(VisitorIndex));
      d2U.resize(d_table->size(VisitorIndex));
      dU.resize(d_table->size(VisitorIndex));

      ValueType sumu = 0.0;
      ValueType dudr, d2udr2;
      for(int i=0; i<d_table->size(SourceIndex); i++) {
	for(int nn=d_table->M[i]; nn<d_table->M[i+1]; nn++) {
	  int j = d_table->J[nn];
	  //U[j] += F[d_table->PairID[nn]]->evaluate(d_table->r(nn));
	  //Grad/Lap are not calculated here
	  U[j] += F[d_table->PairID[nn]]->evaluate(d_table->r(nn), dudr, d2udr2);
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

    /** copy the current data from a buffer
     *@param P the ParticleSet to operate on
     *@param buf PooledData which stores the data for each walker
     *
     *copyFromBuffer uses the data stored by registerData or evaluate(P,buf)
     */
    void copyFromBuffer(ParticleSet& P, PooledData<RealType>& buf) {
      buf.get(U.begin(), U.end());
      buf.get(d2U.begin(), d2U.end());
      buf.get(FirstAddressOfdU,LastAddressOfdU);
    }

    /** return the current value and copy the current data to a buffer
     *@param P the ParticleSet to operate on
     *@param buf PooledData which stores the data for each walker
     */
    inline ValueType evaluate(ParticleSet& P, PooledData<RealType>& buf) {
      ValueType sumu = 0.0;
      for(int i=0; i<U.size(); i++) sumu+=U[i];

      buf.put(U.begin(), U.end());
      buf.put(d2U.begin(), d2U.end());
      buf.put(FirstAddressOfdU,LastAddressOfdU);

      return exp(-sumu);
    }

#ifdef USE_FASTWALKER
    void evaluate(WalkerSetRef& W, 
		  ValueVectorType& psi,
		  WalkerSetRef::WalkerGradient_t& G,
		  WalkerSetRef::WalkerLaplacian_t& L) {

      ValueType dudr, d2udr2;
      int nw = W.walkers();
      const DistanceTableData::IndexVectorType& M = d_table->M;
      const DistanceTableData::IndexVectorType& J = d_table->J;
      const DistanceTableData::IndexVectorType& PairID = d_table->PairID;
      for(int iw=0; iw<nw; iw++) {
	ValueType sumu = 0.0;
	for(int i=0; i<d_table->size(SourceIndex); i++) {
	  for(int nn=M[i]; nn<M[i+1]; nn++) {
	    int j = J[nn];
	    sumu += F[PairID[nn]]->evaluate(d_table->r(iw,nn), dudr, d2udr2);
	    dudr *= d_table->rinv(iw,nn);
	    G(iw,j) -= dudr*d_table->dr(iw,nn);
	    L(iw,j) -= d2udr2+2.0*dudr;
	  }
	}
	psi[iw] *=  exp(-sumu);
      }
    }
#else
    void evaluate(WalkerSetRef& W,
                  ValueVectorType& psi,
                  WalkerSetRef::WalkerGradient_t& G,
                  WalkerSetRef::WalkerLaplacian_t& L) {
      ValueType dudr, d2udr2;
      int nw = W.walkers();
      const DistanceTableData::IndexVectorType& M = d_table->M;
      const DistanceTableData::IndexVectorType& J = d_table->J;
      const DistanceTableData::IndexVectorType& PairID = d_table->PairID;
      vector<ValueType> sumu(nw,0.0);
      for(int i=0; i<d_table->size(SourceIndex); i++) {
        for(int nn=M[i]; nn<M[i+1]; nn++) {
          int j = J[nn];
          for(int iw=0; iw<nw; iw++) {
            sumu[iw] += F[PairID[nn]]->evaluate(d_table->r(iw,nn),dudr,d2udr2);
            dudr *= d_table->rinv(iw,nn);
            G(iw,j) -= dudr*d_table->dr(iw,nn);
            L(iw,j) -= d2udr2+2.0*dudr;
          }
        }
      }
      for(int iw=0; iw<nw; iw++) psi[iw]*= exp(-sumu[iw]);
    }
#endif
  };

}
#endif
/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$ 
 ***************************************************************************/

