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
#ifndef OHMMS_QMC_DISTANCETABLEDATAIMPL_H
#define OHMMS_QMC_DISTANCETABLEDATAIMPL_H

#include "Configuration.h"
#include "Particle/ParticleSet.h"

//// experimental version
namespace ohmmsqmc {

  class ParticleSet;

  /** Class to store pair data between two ParticleSets.
   * @author Jeongnim Kim 
   * @brief DistanceTableData is determined by Source and Target.
   */
  class DistanceTableData {

  public:

    /**enum for index ordering and storage. 
     *@brief Equivalent to using three-dimensional array with (i,j,k)
     * for i = source particle index (slowest), j = target particle index
     * and k = copies (walkers) index.
     */
    enum {SourceIndex=0, VisitorIndex, WalkerIndex};

    ///size of indicies
    TinyVector<int,3> N;

    /**@defgroup 
     * an auxiliary array to handle connections or nearest neighbors
     *@{
     *@brief M.size() = N[SourceIndex]+1
     * M[i+i] - M[i] = the number of connected points to the i-th source
     */ 
    IndexVector_t  M;

    /*@brief J.size() = M[N[SourceIndex]]
     * J[nn] = the index of the connected point for the i-th point 
     * satisfying  \f$M[i] <= nn < M[i+i]\f$
     */ 
    IndexVector_t  J;
    /** @}*/

    /**defgroup storage
     *@{
     *@brief displacement vectors \f$dr(i,j) = R(j)-R(i)\f$ 
     */
    PosVector_t dr;
    
    /*@brief Carteisan distances\f$r(i,j) = |R(j)-R(i)|\f$ */
    ValueVector_t r;
    
    /*@brief Inverse of Carteisan distances\f$rinv(i,j) = 1/|R(j)-R(i)|\f$ */
    ValueVector_t rinv;
    /**@}*/
    
    ///constructor using source and target ParticleSet
    DistanceTableData(const ParticleSet& source, const ParticleSet& target)
    {  }

    ///virutal destructor
    virtual ~DistanceTableData() { }

    ///returns the size of each dimension using enum
    inline int size(int i) const { return N[i];}

    inline int getTotNadj() const { return M[M.size()-1];}

    ///returns the serical index of (i,j,iw)
    inline int index(int i, int j, int iw) const { 
      return N[WalkerIndex]*(j+N[VisitorIndex]*i)+iw;
    }

    ///evaluate the Distance Table using ActiveWalkers
    virtual void 
    evaluate(const PosVector_t& a, int visitors, int copies = 1) = 0;

    ///evaluate the Distance Table using only with position array
    virtual void evaluate(const PosVector_t& a) = 0;

  protected:
    /**resize the storage
     *@param n1 the number of pairs which is evaluated by a derived class
     *@param n2 the number of copies
     *@brief The data for the pair distances, normalized displacements
     *and the distance inverses are stored in a linear storage.
     * The logical view of these storages is (ipair,iwalker),
     * where 0 <= ipair < M[N[SourceIndex]] and 0 <= iwalker < N[WalkerIndex]
     * This scheme can handle both dense and sparse distance tables,
     * and full or half of the pairs.
     * Note that this function is protected and the derived classes are
     * responsible to call this function for memory allocation and any
     * change in the indices N.
     */
    void resize(int n1, int n2=1) {
      N[WalkerIndex] = n2;
      int n = n1*n2;
      dr.resize(n);
      r.resize(n);
      rinv.resize(n);
    }


  private:

    ///disable copy constructor of the base class
    DistanceTableData(const DistanceTableData& a) { }
  };

  /** A derived classe from DistacneTableData, specialized for dense-symmetric case,
   * i.e., the source and target sets are identical.
   *@todo Template with the boundary conditions
   */
  struct SymmetricDTD: public DistanceTableData {
 
    ///constructor using source and target arrays
    SymmetricDTD(const ParticleSet& source, const ParticleSet& target):
      DistanceTableData(source,target){ }

    void reset(int m, int nactive);

    ///evaluate the Distance Table using a set of Particle Positions
    void evaluate(const PosVector_t& a, int visitors, int copies = 1);

    ///not so useful inline but who knows
    inline void evaluate(const PosVector_t& a){
      int n = a.size();
      reset(n,1);
      int ij = 0;
      for(int i=0; i<n; i++) {
	for(int j=i+1; j<n; j++, ij++) {
	  SPPosition_t drij = a[j]-a[i];
	  value_type sep = sqrt(dot(drij,drij));
	  r(ij) = sep;
	  rinv(ij) = 1.0/sep;      
	  dr(ij) = drij;
	}
      }
    }
  };


  /** A derived classe from DistacneTableData, specialized for dense-asymmetric case,
   * i.e., the source and target are distinct.
   *@todo Template with the boundary conditions
   */
  struct AsymmetricDTD: public DistanceTableData {

    const ParticleSet& Origin;
    
    AsymmetricDTD(const ParticleSet& source, 
		  const ParticleSet& target): 
      DistanceTableData(source,target), Origin(source){}
    
    void reset(int n1, int n2, int nactive);

    ///evaluate the Distance Table using a set of Particle Positions
    void evaluate(const PosVector_t& a, int visitors, int copies = 1);

    ///not so useful inline but who knows
    inline void evaluate(const PosVector_t& a){
    
      reset(Origin.getTotalNum(),a.size(),1);
      int ij=0;
      for(int i=0; i<N[SourceIndex]; i++) {
	SPPosition_t r0 = Origin.R[i];
	for(int j=0; j<N[VisitorIndex]; j++,ij++) {
	  SPPosition_t drij = a[j]-r0;
	  value_type sep = sqrt(dot(drij,drij));
	  r(ij)    = sep;
	  rinv(ij) = 1.0/sep;
	  dr(ij)   = drij;
	}
      }
    }
  };

}
#endif
/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$ 
 ***************************************************************************/
