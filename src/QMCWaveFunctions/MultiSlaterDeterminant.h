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
#ifndef OHMMS_QMC_MULTISLATERDETERMINANT_ORBITAL_H
#define OHMMS_QMC_MULTISLATERDETERMINANT_ORBITAL_H
#include "Configuration.h"
#include "QMCWaveFunctions/OrbitalBase.h"
#include "QMCWaveFunctions/SlaterDeterminant.h"

namespace ohmmsqmc {

 /**
   *@brief An AntiSymmetric OrbitalBase composed of a linear
   *combination of SlaterDeterminants. 
   *
   *\f[ 
   *MS({\bf R}) = \sum_n c_n S_n({\bf R}) 
   *\f].
   *
   *The (S)ingle(P)article(O)rbitalSet template parameter is an 
   *engine which fills in single-particle orbital terms.  
   * 
   \f[
   \frac{\nabla_i\Psi}{\Psi} = \frac{\sum_{n=1}^M c_n \nabla_i D_n}
   {\sum_{n=1}^M c_n D_n}
   \f]
   \f[ 
   \frac{\sum_{n=1}^M c_n S_n(\sum_{j=1}^N(\nabla_i
   S^{ij}_n({\bf r_i}))(S^{-1})^{ji}_n}{\sum_{n=1}^M c_n S_n}
   \f]
   The Laplacian
   \f[
   \frac{\nabla^2_i\Psi}{\Psi} = \frac{\sum_{n=1}^M c_n S_n(\sum_{j=1}^N
   (\nabla_i^2S^{ij}_n({\bf r_i}))(S^{-1})^{ji}_n}{\sum_{n=1}^M c_n S_n}
   \f]
   */
template<class SPOSet>
class MultiSlaterDeterminant: public OrbitalBase {

public:

  typedef SlaterDeterminant<SPOSet> DeterminantSet_t;

  ///constructor
  MultiSlaterDeterminant() { Optimizable=false;}

  ///destructor
  ~MultiSlaterDeterminant() { }

  /**
     add a new SlaterDeterminant with coefficient c to the 
     list of determinants
  */
  void add(DeterminantSet_t* sdet, RealType c) {
    SDets.push_back(sdet);
    C.push_back(c);
  }

  void reset() {  
    if(Optimizable) for(int i=0; i<SDets.size(); i++) SDets[i]->reset();
  }

  void initParameters() { }

  inline ValueType
  evaluate(ParticleSet& P, //const DistanceTableData* dtable,
	   ParticleSet::ParticleGradient_t& G,
	   ParticleSet::ParticleLaplacian_t& L){

    int n = P.getTotalNum();
    ParticleSet::ParticleGradient_t g(n), gt(n);
    ParticleSet::ParticleLaplacian_t l(n), lt(n);

    ValueType psi = 0.0;
    for(int i=0; i<SDets.size(); i++){
      ValueType cdet = C[i]*SDets[i]->evaluate(P,g,l);
      psi += cdet;
      gt += cdet*g;
      lt += cdet*l;
      g=0.0;
      l=0.0;
    }
    ValueType psiinv = 1.0/psi;
    G += gt*psiinv;
    L += lt*psiinv;
    return psi;
  }

  inline ValueType
  evaluateLog(ParticleSet& P, //const DistanceTableData* dtable,
	      ParticleSet::ParticleGradient_t& G,
	      ParticleSet::ParticleLaplacian_t& L){
    ValueType psi = evaluate(P,G,L);
    SignValue = (psi<0.0)?-1.0:1.0;
    LogValue = log(abs(psi));
    return LogValue;
  }

  /// returns the dimension of the j-th determinant 
  inline int size(int i, int j) const {return SDets[i]->size(j);}


  inline void evaluate(WalkerSetRef& W, //const DistanceTableData* dtable,
		       ValueVectorType& psi,
		       WalkerSetRef::WalkerGradient_t& G,
		       WalkerSetRef::WalkerLaplacian_t& L) {
    for(int i=0; i<SDets.size(); i++) SDets[i]->evaluate(W,psi,G,L);
      //Dets[i]->evaluate(W,dtable,psi,G,L);
  }
  
  void registerData(ParticleSet& P, PooledData<RealType>& buf){
    std::cerr << "MultiSlaterDeterminant::registerData is empty" << std::endl;
  }
  
  void copyFromBuffer(ParticleSet& P, PooledData<RealType>& buf) {
    std::cerr << "MultiSlaterDeterminant::putData is empty" << std::endl;
  }
  
  inline ValueType
  ratio(ParticleSet& P, int iat, 
	ParticleSet::ParticleGradient_t& G, 
	ParticleSet::ParticleLaplacian_t& L) { 
    std::cerr << "MultiSlaterDeterminant should not be used by Particle-By-Particle update"
              << std::endl;
    return 1.0;
  }

  inline void restore(int iat) { }

  void update(ParticleSet& P, int iat) {
    std::cerr << "MultiSlaterDeterminant::update is empty" << std::endl;
  }

  inline ValueType
  ratio(ParticleSet& P, int iat) { 
    std::cerr << "MultiSlaterDeterminant should not be used by Particle-By-Particle update"
              << std::endl;
    return 1.0;
  }

  void update(ParticleSet& P, 
	      ParticleSet::ParticleGradient_t& G, 
	      ParticleSet::ParticleLaplacian_t& L,
	      int iat) {
    std::cerr << "MultiSlaterDeterminant::update is empty" << std::endl;
  }
  
  ValueType evaluate(ParticleSet& P, PooledData<RealType>& buf) {
    std::cerr << "MultiSlaterDeterminant::evaluate is empty" << std::endl;
    return 1.0;
  }

private:
  vector<DeterminantSet_t*> SDets;
  vector<RealType> C;
};

}
#endif
/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$
 ***************************************************************************/
