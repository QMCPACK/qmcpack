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
#ifndef QMCPLUSPLUS_DISTANCETABLEDATAIMPL_H
#define QMCPLUSPLUS_DISTANCETABLEDATAIMPL_H

#include "Particle/ParticleSet.h"
#include "Utilities/PooledData.h"
#include "OhmmsPETE/OhmmsVector.h"
#include "OhmmsPETE/OhmmsMatrix.h"
#include <bitset>

namespace qmcplusplus {

  /** container for the pair data
   */
  template<class T, unsigned D>
  struct PairDataType {
    ///distance-related data
    T   r, rr, rinv;
    ///displacement vector
    TinyVector<T,D> dr;
    ///default constructor
    inline PairDataType(){}
    ///copy constructor
    inline PairDataType(const PairDataType<T,D>& p):r(p.r),rr(p.rr),rinv(p.rinv),dr(p.dr){}
    ///copy operator
    inline PairDataType<T,D>& operator=(const PairDataType<T,D>& p) {
      r=p.r;rr=p.rr;rinv=p.rinv;dr=p.dr;
      return *this;
    }
    ///set the values
    inline void set(const TinyVector<T,D>& displ, T sep2) {
      r=sqrt(sep2);
      rr=sep2;
      rinv=1.0/r;
      dr=displ;
    }
    ///clear the value
    inline void reset() {
      r=0.0;
    }
  };

  /** @defgroup nnlist Distance-table group
   * @brief class to manage a set of data for distance relations between ParticleSet objects.
   */
  template<class T, unsigned N>
  struct TempDisplacement {
    ///new distance
    T r1;
    ///inverse of the new distance
    T rinv1;
    ///new displacement
    TinyVector<T,N> dr1;
    TinyVector<T,N> dr1_nobox;
    inline TempDisplacement():r1(0.0),rinv1(0.0){}
    inline void reset() {r1=0.0;rinv1=0.0;dr1=0.0;}
    /////old distance
    //T r0;
    /////inverse of old distance
    //T rinv0;
    /////old displacement
    //TinyVector<T,N> dr0;
    //inline TempDisplacement():r1(0.0),rinv1(0.0),r0(0.0),rinv0(0.0) {}
    //inline void reset() {r1=0.0;rinv1=0.0;dr1=0.0;r0=0.0;rinv0=0.0;dr0=0.0;}
  };


  /** @ingroup nnlist
   * @brief Abstract class to manage pair data between two ParticleSets.
   * 
   * Each DistanceTableData object is fined by Source and Target of ParticleSet types.
   */
  class DistanceTableData: public QMCTraits {

  public:
    
    /**enum for index ordering and storage. 
     *@brief Equivalent to using three-dimensional array with (i,j,k)
     * for i = source particle index (slowest), 
     *     j = target particle index
     *     k = copies (walkers) index.
     */
    enum {WalkerIndex=0, SourceIndex, VisitorIndex, PairIndex};

    typedef std::vector<IndexType>       IndexVectorType;
    typedef TempDisplacement<RealType,DIM> TempDistType;
    typedef PooledData<RealType>           BufferType;

    ///true if bound box exists
    bool UseBoundBox;
    ///Index of the particle  with a trial move
    IndexType activePtcl;
    ///size of indicies
    TinyVector<IndexType,4> N;
    ///** Maximum radius */
    //RealType Rmax;
    ///** Maximum square */
    //RealType Rmax2;

    /** @brief M.size() = N[SourceIndex]+1
     *
     * M[i+i] - M[i] = the number of connected points to the i-th source
     */ 
    IndexVectorType M;

    /** @brief J.size() = M[N[SourceIndex]]
     *
     * J[nn] = the index of the connected point for the i-th point 
     * satisfying  \f$M[i] <= nn < M[i+i]\f$
     */ 
    IndexVectorType J;

    /** @brief PairID.size() = M[N[SourceIndex]]
     *
     * PairID[nn] = the index of the connected point for the i-th point 
     * satisfying  \f$PairIDM[i] <= nn < PairID[i+i]\f$
     */ 
    IndexVectorType PairID;

    /** Locator of the pair  */
    std::vector<IndexType> IJ;

    /** @brief A NN relation of all the source particles with respect to an activePtcl
     *
     * This data is for particle-by-particle move.
     * When a MC move is propsed to the activePtcl, the old and new distance relation
     * is stored in Temp. When the move is accepted, the new data replace the old.
     * If the move is rejected, nothing is done and new data will be overwritten.
     */
    std::vector<TempDistType> Temp;
    std::vector<RealType> temp_r;
    std::vector<PosType> temp_dr;

    std::string Name;
    ///constructor using source and target ParticleSet
    DistanceTableData(const ParticleSet& source, const ParticleSet& target)
      : Origin(source), N(0)//, Rmax(1e6), Rmax2(1e12)
    {  }

    ///virutal destructor
    virtual ~DistanceTableData() { }

    ///return the name of table
    inline string getName() const { return Name;}
    ///set the name of table
    inline void setName(const string& tname) { Name = tname;}
    /////set the maximum radius
    //inline void setRmax(RealType rc) { Rmax=rc;Rmax2=rc*rc;}
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
    //virtual void evaluate(const WalkerSetRef& W) = 0;

    ///evaluate the Distance Table using only with position array
    virtual void evaluate(const ParticleSet& P) = 0;

    ///evaluate the temporary pair relations
    virtual void move(const ParticleSet& P, const PosType& rnew, IndexType jat) =0;

    ///evaluate the temporary pair relations
    virtual void moveby(const ParticleSet& P, const PosType& displ, IndexType jat) =0;

    ///evaluate the distance tables with a sphere move
    virtual void moveOnSphere(const ParticleSet& P, const PosType& displ, IndexType jat) =0;

    ///update the distance table by the pair relations
    virtual void update(IndexType jat) = 0;

    ///create storage for nwalkers
    virtual void create(int walkers) = 0;

    /** @brief register nnlist data to buf so that it can copyToBuffer and copyFromBuffer
     *
     * This function is used for particle-by-particle MC methods to register distance-table data
     * to an anonymous buffer.
     */
    inline void registerData(BufferType& buf) {
      RealType* first = &(dr_m[0][0]);
      buf.add(first,first+npairs_m*DIM);
      buf.add(r_m.begin(), r_m.end());
      buf.add(rinv_m.begin(), rinv_m.end());
    }

    /** @brief register nnlist data to buf so that it can copyToBuffer and copyFromBuffer
     *
     * This function is used for particle-by-particle MC methods to register distance-table data
     * to an anonymous buffer.
     */
    inline void updateBuffer(BufferType& buf) {
      RealType* first = &(dr_m[0][0]);
      buf.put(first,first+npairs_m*DIM);
      buf.put(r_m.begin(), r_m.end());
      buf.put(rinv_m.begin(), rinv_m.end());
    }

    /** @brief copy the data to an anonymous buffer
     *
     * Any data that will be used by the next iteration will be copied to a buffer.
     */
    inline void copyToBuffer(BufferType& buf) {
      RealType* first = &(dr_m[0][0]);
      buf.put(first,first+npairs_m*DIM);
      buf.put(r_m.begin(), r_m.end());
      buf.put(rinv_m.begin(), rinv_m.end());
    }

    /** @brief copy the data from an anonymous buffer
     *
     * Any data that is used by the previous iteration will be copied from a buffer.
     */
    inline void copyFromBuffer(BufferType& buf) {
      RealType* first = &(dr_m[0][0]);
      buf.get(first,first+npairs_m*DIM);
      buf.get(r_m.begin(), r_m.end());
      buf.get(rinv_m.begin(), rinv_m.end());
    }

    inline void print(std::ostream& os) {
      os << "Table " << Origin.getName() << endl; 
      for(int i=0; i<r_m.size(); i++)
	os << r_m[i] << " ";
      os << endl;
    }

  protected:

    const ParticleSet& Origin;

    ///number of pairs
    int npairs_m;

    /**defgroup storage data for nearest-neighbor relations
     */
    /*@{*/
    /** Cartesian distance \f$r(i,j) = |R(j)-R(i)|\f$ */
    std::vector<RealType> r_m;
    /** Cartesian distance \f$rinv(i,j) = 1/r(i,j)\f$ */
    std::vector<RealType> rinv_m;
    /** displacement vectors \f$dr(i,j) = R(j)-R(i)\f$  */
    std::vector<PosType> dr_m;
    /*@}*/

    Matrix<PosType> dr2_m;
    Matrix<RealType> r2_m, rinv2_m;

    /**resize the storage
     *@param npairs number of pairs which is evaluated by a derived class
     *@param nw number of copies
     *
     * The data for the pair distances, displacements
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
	//rr_m.resize(npairs);
	rinv_m.resize(npairs);

	Temp.resize(N[SourceIndex]);
        temp_r.resize(N[SourceIndex]);
        temp_dr.resize(N[SourceIndex]);
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

  };
}
#endif
/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$ 
 ***************************************************************************/
