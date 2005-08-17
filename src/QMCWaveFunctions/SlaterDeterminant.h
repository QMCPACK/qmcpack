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
#ifndef OHMMS_QMC_SLATERDETERMINANT_ORBITAL_H
#define OHMMS_QMC_SLATERDETERMINANT_ORBITAL_H
#include "Configuration.h"
#include "QMCWaveFunctions/OrbitalBase.h"
#include "QMCWaveFunctions/DiracDeterminant.h"

namespace ohmmsqmc {

  /** @ingroup OrbitalComponent
   *  @brief An AntiSymmetric OrbitalBase composed of DiracDeterminants. 
   *
   A SlaterDeterminant is a product of DiracDeterminants
   \f[ 
   S({\bf R}) = \prod_i D_i({\bf r_{first,i},\ldots,r_{last,i}}).
   \f]
   *Typically \f$ S({\bf R}) \f$ is a product of an up and down 
   *DiracDeterminant such that
   \f[
   S({\bf R}) = D_{\uparrow}({\bf r_1,r_2,\ldots,r_{N_{\uparrow}}})
   D_{\downarrow}({\bf r_{N_{\uparrow +1}},\ldots,r_{N}})
   \f]
   *
   *SlaterDeterminant is a composite class which collects each of 
   *the compotents in a simple manner: multiplications of determinants
   *and addition of gradients and laplacian terms in a linear pattern.
   *
   *The (S)ingle(P)article(O)rbitalSet template parameter is an 
   *engine which fills in single-particle orbital terms.  
   * 
   *@note MultiSlaterDeterminant is a linear combination of SlaterDeterminants.
   */
  template<class SPOSet>
  class SlaterDeterminant: public OrbitalBase {

  public:

    typedef DiracDeterminant<SPOSet> Determinant_t;

    /// constructor
    SlaterDeterminant() {M.resize(3,0);Optimizable=false;}

    ///destructor
    ~SlaterDeterminant() { }

    ///add a new DiracDeterminant to the list of determinants
    void add(Determinant_t* det) { 
      int last=Dets.size();
      Dets.push_back(det);
      M[last+1]=M[last]+Dets[last]->rows();
    }

    ///reset all the Dirac determinants, Optimizable is true
    void reset() {  
      if(Optimizable) for(int i=0; i<Dets.size(); i++) Dets[i]->reset();
    }

    /** Calculate the value of the Slater determinant for the input configuration. 
     *@param P input configuration containing N particles
     *@param G a vector containing N gradients
     *@param L a vector containing N laplacians
     *@return SlaterDeterminant value
     *
     *Add the gradient and laplacian contribution of the Slater determinant to G(radient) and L(aplacian)
     *for local energy calculations.
     */
    inline ValueType 
    evaluate(ParticleSet& P, 
	     ParticleSet::ParticleGradient_t& G, 
	     ParticleSet::ParticleLaplacian_t& L) {
      ValueType psi = 1.0;
      for(int i=0; i<Dets.size(); i++) psi *= Dets[i]->evaluate(P,G,L);
      return psi;
    }

    inline ValueType 
    evaluateLog(ParticleSet& P, 
	        ParticleSet::ParticleGradient_t& G, 
	        ParticleSet::ParticleLaplacian_t& L) {
      ValueType psi = 1.0;
      for(int i=0; i<Dets.size(); i++) psi *= Dets[i]->evaluate(P,G,L);
      SignValue = (psi<0.0)?-1.0:1.0;
      LogValue = log(abs(psi));
      return LogValue;
    }

    virtual void resizeByWalkers(int nwalkers) {
      for(int i=0; i<Dets.size(); i++) 	Dets[i]->resizeByWalkers(nwalkers);
    }

    ///return the total number of Dirac determinants
    inline int size() const { return Dets.size();}

    ///return the dimension of the i-th Dirac determinant
    inline int size(int i) const { return Dets[i]->cols();}

    inline void evaluate(WalkerSetRef& W, //const DistanceTableData* dtable, 
			 ValueVectorType& psi,
			 WalkerSetRef::WalkerGradient_t& G,
			 WalkerSetRef::WalkerLaplacian_t& L) {
      for(int i=0; i<Dets.size(); i++) 	Dets[i]->evaluate(W,psi,G,L);
    }

    /** similar to evaluateLog 
     */
    ValueType registerData(ParticleSet& P, PooledData<RealType>& buf){
      //for(int i=0; i<Dets.size(); i++) Dets[i]->registerData(P,buf);
      ValueType psi = 1.0;
      for(int i=0; i<Dets.size(); i++) 
        psi *= Dets[i]->registerData(P,buf);
      SignValue = (psi<0.0)?-1.0:1.0;
      LogValue = log(abs(psi));
      return LogValue;
    }
    
    void copyFromBuffer(ParticleSet& P, PooledData<RealType>& buf) {
      for(int i=0; i<Dets.size(); i++) 	Dets[i]->copyFromBuffer(P,buf);
    }

    /** reimplements the virtual function
     *
     * The DiractDeterminants of SlaterDeterminant need to save the inverse
     * of the determinant matrix to evaluate ratio
     */
    void dumpToBuffer(ParticleSet& P, PooledData<RealType>& buf) {
      for(int i=0; i<Dets.size(); i++) 	Dets[i]->dumpToBuffer(P,buf);
    }

    /** reimplements the virtual function
     *
     * Matching function to dumpToBuffer.
     */
    void dumpFromBuffer(ParticleSet& P, PooledData<RealType>& buf) {
      for(int i=0; i<Dets.size(); i++) 	Dets[i]->dumpFromBuffer(P,buf);
    }

    ValueType evaluate(ParticleSet& P, PooledData<RealType>& buf) {
      ValueType r=1.0;
      for(int i=0; i<Dets.size(); i++) 	r *= Dets[i]->evaluate(P,buf);
      return r;
    }

    inline ValueType ratio(ParticleSet& P, int iat,
			   ParticleSet::ParticleGradient_t& dG, 
			   ParticleSet::ParticleLaplacian_t& dL) { 
      int i=1;
      while(iat>=M[i]) {i++;}
      return  Dets[i-1]->ratio(P,iat, dG, dL);
    }
    
    inline void restore(int iat) {
      int i=1;
      while(iat>=M[i]) {i++;}
      Dets[i-1]->restore(iat);
    }

    inline void update(ParticleSet& P, int iat) {
      int i=1;
      while(iat>=M[i]) {i++;}
      Dets[i-1]->update(P,iat);
    }

    ValueType
    ratio(ParticleSet& P, int iat) {
      int i=1;
      while(iat>=M[i]) {i++;}
      return Dets[i-1]->ratio(P,iat);
    } 	  

    void update(ParticleSet& P, 
		ParticleSet::ParticleGradient_t& dG, 
		ParticleSet::ParticleLaplacian_t& dL,
		int iat) {
      int i=1;
      while(iat>=M[i]) {i++;}
      Dets[i-1]->update(P,dG,dL,iat);
    }


  private:
    vector<int> M;
    ///container for the DiracDeterminants
    vector<Determinant_t*>  Dets;
  };
}
#endif
/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$ 
 ***************************************************************************/
