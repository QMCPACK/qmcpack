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

#include "Particle/ParticleSet.h"
#include "Particle/WalkerSetRef.h"
#include "Utilities/PooledData.h"
#include "OhmmsPETE/OhmmsMatrix.h"

namespace ohmmsqmc {

  /** Class to store pair data between two ParticleSets.
   * @author Jeongnim Kim 
   * @brief DistanceTableData is determined by Source and Target.
   */
  class DistanceTableData: public QMCTraits {

  public:

    
    /**enum for index ordering and storage. 
     *@brief Equivalent to using three-dimensional array with (i,j,k)
     * for i = source particle index (slowest), 
     *     j = target particle index
     *     k = copies (walkers) index.
     */
    enum {WalkerIndex=0, SourceIndex, VisitorIndex};

    typedef std::vector<IndexType>  IndexVectorType;

    ///size of indicies
    TinyVector<IndexType,3> N;

    /**@defgroup 
     * an auxiliary array to handle connections or nearest neighbors
     *@{
     *@brief M.size() = N[SourceIndex]+1
     * M[i+i] - M[i] = the number of connected points to the i-th source
     */ 
    IndexVectorType M;

    /*@brief J.size() = M[N[SourceIndex]]
     * J[nn] = the index of the connected point for the i-th point 
     * satisfying  \f$M[i] <= nn < M[i+i]\f$
     */ 
    IndexVectorType J;


    /*@brief PairID.size() = M[N[SourceIndex]]
     * PairID[nn] = the index of the connected point for the i-th point 
     * satisfying  \f$PairIDM[i] <= nn < PairID[i+i]\f$
     */ 
    IndexVectorType PairID;
    /** @}*/

    ///constructor using source and target ParticleSet
    DistanceTableData(const ParticleSet& source, const ParticleSet& target)
      :Origin(source)
    {  }

    ///virutal destructor
    virtual ~DistanceTableData() { }

    ///returns the reference the origin
    const ParticleSet& origin() const { return Origin;}

    //@{access functions to the distance, inverse of the distance and directional consine vector
#ifdef USE_FASTWALKER
    inline PosType dr(int iw, int iat) const {return dr2_m(iat,iw);}
    inline RealType r(int iw, int iat) const {return r2_m(iat,iw);}
    inline RealType rinv(int iw, int iat) const {return rinv2_m(iat,iw);}
#else
    inline PosType dr(int iw, int iat) const {return dr2_m(iw,iat);}
    inline RealType r(int iw, int iat) const {return r2_m(iw,iat);}
    inline RealType rinv(int iw, int iat) const {return rinv2_m(iw,iat);}
#endif
    inline PosType dr(int j) const {return dr_m[j];}
    inline RealType r(int j) const {return r_m[j];}
    inline RealType rinv(int j) const {return rinv_m[j];}
    //@}

    ///returns the number of centers
    inline IndexType centers() const { return Origin.getTotalNum();}

    ///returns the size of each dimension using enum
    inline IndexType size(int i) const { return N[i];}

    //inline IndexType getTotNadj() const { return M[M.size()-1];}
    inline IndexType getTotNadj() const { return npairs_m;}

    //!< Returns a number of neighbors of the i-th ptcl.
    inline IndexType nadj(int i) const { return M[i+1]-M[i];}
    
    //!< Returns the id of j-th neighbor for i-th ptcl
    inline IndexType iadj(int i, int j) const { return J[M[i] +j];}

    //!< Returns the id of j-th neighbor for i-th ptcl
    inline IndexType loc(int i, int j) const { return M[i] + j;}

    ///evaluate the Distance Table using ActiveWalkers
    virtual void evaluate(const WalkerSetRef& W) = 0;

    ///evaluate the Distance Table using only with position array
    virtual void evaluate(const ParticleSet& P) = 0;

    virtual void update(const ParticleSet& P, IndexType iat) = 0;

    virtual void create(int walkers) = 0;

    //     virtual void select(int iat) =0;
    //     virtual void update(int iat) =0;

    inline void registerData(PooledData<RealType>& buf) {
      RealType* first = &(dr_m[0][0]);
      buf.add(first,first+npairs_m*DIM);
      buf.add(r_m.begin(), r_m.end());
      buf.add(rinv_m.begin(), rinv_m.end());
    }

    inline void getData(PooledData<RealType>& buf) {
      RealType* first = &(dr_m[0][0]);
      buf.get(first,first+npairs_m*DIM);
      buf.get(r_m.begin(), r_m.end());
      buf.get(rinv_m.begin(), rinv_m.end());
    }

    inline void putData(PooledData<RealType>& buf) {
      RealType* first = &(dr_m[0][0]);
      buf.put(first,first+npairs_m*DIM);
      buf.put(r_m.begin(), r_m.end());
      buf.put(rinv_m.begin(), rinv_m.end());
    }

  protected:

    const ParticleSet& Origin;

    ///number of pairs
    int npairs_m;

    /**defgroup storage
     *@{
     *@brief displacement vectors \f$dr(i,j) = R(j)-R(i)\f$ 
     *Inverse of Carteisan distances\f$rinv(i,j) = 1/|R(j)-R(i)|\f$
     */
    std::vector<PosType> dr_m;
    std::vector<RealType> r_m, rinv_m;
    Matrix<PosType> dr2_m;
    Matrix<RealType> r2_m, rinv2_m;
    /*@}*/

    /**resize the storage
     *@param npairs the number of pairs which is evaluated by a derived class
     *@param nw the number of copies
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
    void resize(int npairs, int nw=1) {
      N[WalkerIndex] =nw;
      if(nw==1) {
	dr_m.resize(npairs);
	r_m.resize(npairs);
	rinv_m.resize(npairs);
      } else {
#ifdef USE_FASTWALKER
	dr2_m.resize(npairs,nw);
	r2_m.resize(npairs,nw);
	rinv2_m.resize(npairs,nw);
#else
	dr2_m.resize(nw,npairs);
	r2_m.resize(nw,npairs);
	rinv2_m.resize(nw,npairs);
#endif
      }
    }


  private:
    
    ///disable copy constructor of the base class
    DistanceTableData(const DistanceTableData& a):Origin(a.Origin) { }

  };
}
#endif
/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$ 
 ***************************************************************************/
