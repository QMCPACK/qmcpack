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
#ifndef OHMMS_QMC_LINEARCOMIBINATIONORBITALS_H
#define OHMMS_QMC_LINEARCOMIBINATIONORBITALS_H

#include "OhmmsPETE/OhmmsMatrix.h"
#include "Numerics/DeterminantOperators.h"
#include <limits>

namespace ohmmsqmc {

  class DistanceTableData;

  /** class to handle linear combinations of basis orbitals used to evaluate the Dirac determinants.
   *
   *LCOrbitals stands for (L)inear(C)ombinationOrbitals
   *Any single-particle orbital \f$ \psi_j \f$ that can be represented by
   \f[
   \psi_j ({\bf r}_i) = \sum_k C_{jk} \phi_{k}({\bf r}_i),
   \f]
   *where \f$ \phi_{k} \f$ is the kth basis.
   *The initialization and evaluation of \f$ \phi\ f$'s
   *are delagated to the template class BS, but the linear-combination, i.e.,
   *matrix-vector or matrix-matrix multiplications are handedled by LCOrbitals
   *
   *An example of template class (B)asis(S)et is MolecularOrbitalBasis. 
   */
  template<class BS>
  class LCOrbitals: public OhmmsElementBase {

  public:

    typedef typename BS::RealType  RealType;
    typedef typename BS::ValueType ValueType;
    typedef typename BS::PosType   PosType;
    typedef typename BS::GradType  GradType;

    ///true if the coefficient matrix is the identity matrix
    bool Identity;
    ///id of this object, if ID == 0, the object is responsble for the BasisSet
    int ID;
    ///number of particles handled by this object
    int NumPtcls;
    ///size of basis set
    int BasisSize;
    ///pointer to the basis set
    BS* BasisSet;
    ///matrix containing the coefficients
    Matrix<ValueType> C;

    inline LCOrbitals(BS* bs, int id): 
      Identity(false),ID(id), NumPtcls(0), BasisSet(bs){ }

    inline ~LCOrbitals() {if(!ID) delete BasisSet;}

    ///set the distance table for all the basis sets
    inline void setTable(DistanceTableData* atable) { 
      BasisSet->setTable(atable);
    }

    ///return the number of single particle orbitals
    inline int numOrbitals() const { return C.rows();}

    ///return the number of basis functions
    inline int numBasis() const { return BasisSet->TotalBasis;}

    inline void reset() { BasisSet->reset();}

    ///resize the internal storage of BasisSet by the number of particles
    inline void resize(int nptcl) {
      NumPtcls = nptcl;
      BasisSet->resize(nptcl);
      BasisSize=BasisSet->TotalBasis;
    }

    /** resize the internal storage of BasisSet by the number of walkers
     *
     *@todo This works only with MolecularOrbitalBasis at this point.
     *To be generalized for any basis function
     */
    inline void resizeByWalkers(int nw) {
      BasisSet->resizeByWalkers(nw);
    }

    /** evaluate the values of the single-particle orbitals 
     *@param P input configuration containing N particles
     *@param iat particle index whose position has been modified
     *@param jorb psi values of the single-particle orbitals
     *
     * This function is introduced to function evaluations.
     */
    inline ValueType
    evaluate(const ParticleSet& P, int iat, int jorb) {
      BasisSet->evaluate(P,iat);
      return dot(&C(jorb,0),BasisSet->y(0),BasisSize);
    }

    /** evaluate the values of the single-particle orbitals 
     *@param P input configuration containing N particles
     *@param iat particle index
     *@param psi values of the single-particle orbitals
     *
     * This function is introduced to function evaluations.
     */
    template<class VV>
    inline void 
    evaluate(const ParticleSet& P, int iat, VV& psi) {
      BasisSet->evaluate(P,iat);
      for(int j=0 ; j<NumPtcls; j++) 
        psi[j] = dot(&C(j,0),BasisSet->y(0),BasisSize);
    }

    /**@ingroup particlebyparticle 
     *@brief evaluate the values of the single-particle orbitals for the iat-th particle
     *@param P input configuration containing N particles
     *@param iat particle index
     *@param psi values of the single-particle orbitals
     *@param dpsi gradients of the single-particle orbitals
     *@param d2psi laplacians of the single-particle orbitals
     *
     * This function completes a row of a Dirac Deterimant to perform ratio/update.
     * The particle index identifies the particle whose position is updated.
     */
    template<class VV, class GV>
    inline void 
    evaluate(const ParticleSet& P, int iat, VV& psi, GV& dpsi, VV& d2psi) {
      BasisSet->evaluateAll(P,iat);
      if(Identity) {
        for(int j=0 ; j<NumPtcls; j++) {
	  psi[j] = BasisSet->Y(0,j);
	  dpsi[j] = BasisSet->dY(0,j);
          d2psi[j] = BasisSet->d2Y(0,j);
        }
      } else {
	//int nb = BasisSet->TotalBasis;
	for(int j=0 ; j<NumPtcls; j++) {
	  psi[j]   = dot(&C(j,0),BasisSet->y(0),BasisSize);
	  dpsi[j]  = dot(&C(j,0),BasisSet->dy(0),BasisSize);
	  d2psi[j] = dot(&C(j,0),BasisSet->d2y(0),BasisSize);
        }
      }
    }

    /** complete the values of the single-particle orbitals and their gradients and laplacians
     *@param P input configuration containing N particles
     *@param first index of the first particle
     *@param last index of the last particle
     *@param logdet \f$ logdet(j,i) = \sum_k C_{jk} \phi_{k}({\bf r}_i) \f$
     *@param dlogdet \f$ dlogdet(i,j) =\sum_k C_{jk} \nabla_i \phi_{k}({\bf r}_i) \f$
     *@param d2logdet \f$ d2logdet(i,j) = \sum_k C_{jk} \nabla^2_i \phi_{k}({\bf r}_i) \f$
     *
     *Each object handles [first,last) quantum particles. 
     *The physical particle index maps to the orbital index of [0,last-first).
     */
    template<class VM, class GM>
    inline void 
    evaluate(const ParticleSet& P, int first, int last,
	     VM& logdet, GM& dlogdet, VM& d2logdet) {
      if(!(ID || first)) BasisSet->evaluate(P);
      NumPtcls=last-first;
      //check if identity matrix
      if(Identity) {
	for(int i=0; i<NumPtcls; i++,first++){
	  for(int j=0 ; j<NumPtcls; j++) {
	    logdet(j,i) = BasisSet->Y(first,j);
	    dlogdet(i,j) = BasisSet->dY(first,j);
	    d2logdet(i,j) = BasisSet->d2Y(first,j);
	  }
	}
      } else {
	//int nb = BasisSet->TotalBasis;
	//iat is an index offset for the particle number
	for(int i=0, iat=first; i<NumPtcls; i++,iat++){
	  for(int j=0 ; j<NumPtcls; j++) {
	    //logdet(j,i) = \f$\sum_k^{nb} C(j,k) Y(i+first,k)\f$
	    logdet(j,i) = dot(&C(j,0),BasisSet->y(iat),BasisSize);
	    dlogdet(i,j) = dot(&C(j,0),BasisSet->dy(iat),BasisSize);
	    d2logdet(i,j) = dot(&C(j,0),BasisSet->d2y(iat),BasisSize);
	  }
	}
      }
    }

    template<class VM, class GM>
    inline void 
    evaluate(const WalkerSetRef& W, int first, int last,
	     vector<VM>& logdet, vector<GM>& dlogdet, vector<VM>& d2logdet) {
      //calculate everything
      if(!(ID || first)) {
	BasisSet->evaluate(W);
      }

      int nptcl = last-first;
      if(Identity) {
	for(int iw=0; iw<W.walkers(); iw++) {
	  for(int i=0, iat=first; i<nptcl; i++,iat++){
	    for(int j=0 ; j<nptcl; j++) {
	      logdet[iw](j,i) = *(BasisSet->y(iw,iat));
	      dlogdet[iw](i,j) = *(BasisSet->dy(iw,iat));
	      d2logdet[iw](i,j) = *(BasisSet->d2y(iw,iat));
	    }
	  }
	}
      } else {
	//int nb = BasisSet->TotalBasis;
	int ntot = W.particles();
        for(int iw=0; iw<W.walkers(); iw++) {
  	  for(int i=0, iat=first; i<nptcl; i++, iat++){
	    for(int j=0 ; j<nptcl; j++) {
	      logdet[iw](j,i) = dot(&C(j,0),BasisSet->y(iw,iat),BasisSize);
	      dlogdet[iw](i,j) = dot(&C(j,0),BasisSet->dy(iw,iat),BasisSize);
	      d2logdet[iw](i,j) = dot(&C(j,0),BasisSet->d2y(iw,iat),BasisSize);
	    }
	  }
	}
      }
    }

    /** Parse the xml file for information on the Dirac determinants.
     *@param cur the current xmlNode
     */
    bool put(xmlNodePtr cur) {

      int norb=atoi((const char*)(xmlGetProp(cur, (const xmlChar *)"orbitals")));

      vector<RealType> occupation(norb,1.0);
      vector<ValueType> Ctemp;

      int nocc(0),total(norb);
      cur = cur->xmlChildrenNode;
      Identity = true;

      XMLReport("The number of orbitals for a Dirac Determinant " << norb)
      XMLReport("The number of basis functions " << numBasis())
      while(cur != NULL) {
	string cname((const char*)(cur->name));
	if(cname == "occupation") {
	  nocc = atoi((const char*)(xmlGetProp(cur, (const xmlChar *)"size")));
	  occupation.resize(total);	
	  putContent(occupation,cur);
	} else if(cname == "coefficient" || cname == "parameter" || cname == "Var") {
	  if(xmlHasProp(cur, (const xmlChar*)"size")){
	    total = atoi((const char*)(xmlGetProp(cur, (const xmlChar *)"size")));
	  } 
	  Ctemp.resize(total*numBasis());
	  putContent(Ctemp,cur);
	  Identity = false;
	}
	cur = cur->next;
      }

      if(nocc && nocc != total) {
        ERRORMSG("Inconsistent input for the occupation and size of the coefficients")
        Identity=true;
      }
      cout.setf(ios::scientific, ios::floatfield);
      cout.setf(ios::right,ios::adjustfield);
      cout.precision(12);

      C.resize(norb,numBasis());
      if(Identity) {
        C=0.0;
        for(int i=0; i<norb; i++) C(i,i)=1.0;
      } else {
	int n=0,i=0;
        typename vector<ValueType>::iterator cit(Ctemp.begin());
        BasisSize=numBasis();
	while(i<norb){
	  if(occupation[n]>numeric_limits<RealType>::epsilon()){
            std::copy(cit,cit+BasisSize,C[i]);
	    i++; 
	  }
	  n++;cit+=BasisSize;
	}
      }

      //XMLReport("The Molecular Orbital Coefficients:")
      //XMLReport(C)
      return true;
    }

    ///write to a ostream
    bool get(std::ostream& ) const { return true;}

    ///read from istream
    bool put(std::istream& ) { return true;}

  };
}
#endif

/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$ 
 ***************************************************************************/

