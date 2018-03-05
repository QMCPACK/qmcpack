//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Jaron T. Krogel, krogeljt@ornl.gov, Oak Ridge National Laboratory
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////
    
    



#ifndef QMCPLUSPLUS_DISTANCETABLEDATAIMPL_H
#define QMCPLUSPLUS_DISTANCETABLEDATAIMPL_H

#include "Particle/ParticleSet.h"
#include "Utilities/PooledData.h"
#include "OhmmsPETE/OhmmsVector.h"
#include "OhmmsPETE/OhmmsMatrix.h"
#include "simd/allocator.hpp"
#include <OhmmsSoA/VectorSoaContainer.h>
#include <limits>
#include <bitset>

namespace qmcplusplus
{

/** container for the pair data
 */
template<class T, unsigned D>
struct PairDataType
{
  ///distance-related data
  T   r, rr, rinv;
  ///displacement vector
  TinyVector<T,D> dr;
  ///default constructor
  inline PairDataType() {}
  ///copy constructor
  inline PairDataType(const PairDataType<T,D>& p):r(p.r),rr(p.rr),rinv(p.rinv),dr(p.dr) {}
  ///copy operator
  inline PairDataType<T,D>& operator=(const PairDataType<T,D>& p)
  {
    r=p.r;
    rr=p.rr;
    rinv=p.rinv;
    dr=p.dr;
    return *this;
  }
  ///set the values
  inline void set(const TinyVector<T,D>& displ, T sep2)
  {
    r=sqrt(sep2);
    rr=sep2;
    rinv=1.0/r;
    dr=displ;
  }
  ///clear the value
  inline void reset()
  {
    r=0.0;
  }
};

/** @defgroup nnlist Distance-table group
 * @brief class to manage a set of data for distance relations between ParticleSet objects.
 */
template<class T, unsigned N>
struct TempDisplacement
{
  ///new distance
  T r1;
  ///inverse of the new distance
  T rinv1;
  ///new displacement
  TinyVector<T,N> dr1;
  inline TempDisplacement():r1(0.0),rinv1(0.0) {}
  inline void reset()
  {
    r1=0.0;
    rinv1=0.0;
    dr1=0.0;
  }
};

/** enumerator for DistanceTableData::DTType
 *
 * - DT_AOS Use original AoS type
 * - DT_SOA Use SoA type
 * - DT_AOS_PREFERRED Create AoS type, if possible.
 * - DT_SOA_PREFERRED Create SoA type, if possible.
 * The first user of each pair will decide the type of distance table.
 * It is the responsibility of the user class to check DTType.
 */
enum DistTableType {DT_AOS=0,  DT_SOA, DT_AOS_PREFERRED, DT_SOA_PREFERRED};

/** @ingroup nnlist
 * @brief Abstract class to manage pair data between two ParticleSets.
 *
 * Each DistanceTableData object is fined by Source and Target of ParticleSet types.
 */
struct DistanceTableData 
{
  CONSTEXPR static unsigned DIM=OHMMS_DIM;

  /**enum for index ordering and storage.
   *@brief Equivalent to using three-dimensional array with (i,j,k)
   * for i = source particle index (slowest),
   *     j = target particle index
   *     k = copies (walkers) index.
   */
  enum {WalkerIndex=0, SourceIndex, VisitorIndex, PairIndex};
#if (__cplusplus >= 201103L)
  using IndexType=QMCTraits::IndexType;
  using RealType=QMCTraits::RealType;
  using PosType=QMCTraits::PosType;
  using IndexVectorType=aligned_vector<IndexType>;
  using TempDistType=TempDisplacement<RealType,DIM>;
  using ripair=std::pair<RealType,IndexType>;
  using RowContainer=VectorSoaContainer<RealType,DIM>;
#else
  typedef QMCTraits::IndexType IndexType;
  typedef QMCTraits::RealType RealType;
  typedef QMCTraits::PosType PosType;
  typedef aligned_vector<IndexType> IndexVectorType;
  typedef TempDisplacement<RealType,DIM> TempDistType;
  typedef std::pair<RealType,IndexType> ripair;
  typedef VectorSoaContainer<RealType,DIM> RowContainer;
#endif

  ///type of cell
  int CellType;
  ///ID of this table among many
  int ID;
  ///Type of DT
  int DTType;
  ///size of indicies
  TinyVector<IndexType,4> N;
  ///true, if ratio computations need displacement, e.g. LCAO type
  bool NeedDisplacement;
  ///Maximum radius set by a Hamiltonian
  RealType Rmax;
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
  IndexVectorType IJ;

  /** @brief A NN relation of all the source particles with respect to an activePtcl
   *
   * This data is for particle-by-particle move.
   * When a MC move is propsed to the activePtcl, the old and new distance relation
   * is stored in Temp. When the move is accepted, the new data replace the old.
   * If the move is rejected, nothing is done and new data will be overwritten.
   */
  std::vector<TempDistType> Temp;

  /**defgroup SoA data */
  /*@{*/
  /** Distances[i][j] , [Nsources][Ntargets] */
  Matrix<RealType, aligned_allocator<RealType> > Distances;

  /** Displacements[Nsources]x[3][Ntargets] */
  std::vector<RowContainer> Displacements;

  ///actual memory for Displacements
  aligned_vector<RealType> memoryPool;

  /** temp_r */
  aligned_vector<RealType> Temp_r;

  /** temp_dr */
  RowContainer Temp_dr;

  /** true, if full table is needed at loadWalker */
  bool Need_full_table_loadWalker;
  /*@}*/

  ///name of the table
  std::string Name;
  ///constructor using source and target ParticleSet
  DistanceTableData(const ParticleSet& source, const ParticleSet& target)
    : Origin(&source), N(0), NeedDisplacement(false), Need_full_table_loadWalker(false)
  {  
    Rmax=0; //set 0
    //Rmax=source.Lattice.WignerSeitzRadius;   
  }

  ///virutal destructor
  virtual ~DistanceTableData() { }

  ///return the name of table
  inline std::string getName() const
  {
    return Name;
  }
  ///set the name of table
  inline void setName(const std::string& tname)
  {
    Name = tname;
  }
  ///set the maximum radius
  inline void setRmax(RealType rc) { Rmax=std::max(Rmax,rc);}

  ///returns the reference the origin particleset
  const ParticleSet& origin() const
  {
    return *Origin;
  }
  inline void reset(const ParticleSet* newcenter)
  {
    Origin=newcenter;
  }

  inline bool is_same_type(int dt_type) const
  {
    return DTType == dt_type;
  }

  //@{access functions to the distance, inverse of the distance and directional consine vector
  inline PosType dr(int j) const
  {
    return dr_m[j];
  }
  inline RealType r(int j) const
  {
    return r_m[j];
  }
  inline RealType rinv(int j) const
  {
    return rinv_m[j];
  }

  //@}

  ///returns the number of centers
  inline IndexType centers() const
  {
    return Origin->getTotalNum();
  }
      
  ///returns the number of centers
  inline IndexType targets() const
  {
    return N[VisitorIndex];
  }
  
  ///returns the size of each dimension using enum
  inline IndexType size(int i) const
  {
    return N[i];
  }

  inline IndexType getTotNadj() const
  {
    return npairs_m;
  }

  /// return the distance |R[iadj(i,nj)]-R[i]|
  inline RealType distance(int i, int nj) const
  {
    return (DTType)? r_m2(i,nj): r_m[M[i]+nj];
  }

  /// return the displacement R[iadj(i,nj)]-R[i]
  inline PosType displacement(int i, int nj) const
  {
    return (DTType)? dr_m2(i,nj): dr_m[M[i]+nj];
  }

  //!< Returns a number of neighbors of the i-th ptcl.
  inline IndexType nadj(int i) const
  {
    return (DTType)? M[i]:M[i+1]-M[i];
  }
  //!< Returns the id of nj-th neighbor for i-th ptcl
  inline IndexType iadj(int i, int nj) const
  {
    return (DTType)? J2(i,nj):J[M[i] +nj];
  }
  //!< Returns the id of j-th neighbor for i-th ptcl
  inline IndexType loc(int i, int j) const
  {
    return M[i] + j;
  }

  /** search the closest source particle within rcut
   * @param rcut cutoff radius
   * @return the index of the closest particle
   *
   * Return -1 if none is within rcut
   * @note It searches the temporary list after a particle move is made
   */
  inline IndexType find_closest_source(RealType rcut) const
  {
    int i=0;
    while(i<N[SourceIndex])
    {
      if(Temp[i].r1<rcut) return i;
      i++;
    }
    return -1;
  }

  /** search the closest source particle of the iel-th particle within rcut 
   * @param rcut cutoff radius
   * @return the index of the first particle
   *
   * Return -1 if none is within rcut
   * @note Check the real distance table and only works for the AsymmetricDistanceTable
   */
  inline IndexType find_closest_source(int iel, RealType rcut) const
  {
    int i=0,nn=iel;
    while(nn<r_m.size())
    {
      if(r_m[nn]<rcut) return i;
      nn+=N[VisitorIndex];
      i++;
    }
    return -1;
  }

  ///prepare particle-by-particle moves
  virtual void setPbyP() { }

  ///update internal data after completing particle-by-particle moves
  virtual void donePbyP() { }

  ///evaluate the Distance Table using only with position array
  virtual void evaluate(ParticleSet& P) = 0;

  /// evaluate the Distance Table
  virtual void evaluate(ParticleSet& P, int jat)=0;

  ///evaluate the temporary pair relations
  virtual void move(const ParticleSet& P, const PosType& rnew) =0;

  ///evaluate the distance tables with a sphere move
  virtual void moveOnSphere(const ParticleSet& P, const PosType& rnew) =0;

  ///update the distance table by the pair relations
  virtual void update(IndexType jat) = 0;

  /** build a compact list of a neighbor for the iat source
   * @param iat source particle id
   * @param rcut cutoff radius
   * @param jid compressed index
   * @param dist compressed distance
   * @param displ compressed displacement
   * @return number of target particles within rcut
   */
  virtual size_t get_neighbors(int iat, RealType rcut, int* restrict jid, RealType* restrict dist, PosType* restrict displ) const
  {
    return 0;
  }

  /** find the nearest neighbor
   * @param iat source particle id
   * @return the id of the nearest particle, -1 not found
   */
  virtual int get_first_neighbor(IndexType iat, RealType& r, PosType& dr, bool newpos) const
  {
    return -1;
  }

  /** build a compact list of a neighbor for the iat source
   * @param iat source particle id
   * @param rcut cutoff radius
   * @param dist compressed distance
   * @return number of target particles within rcut
   */
  virtual size_t get_neighbors(int iat, RealType rcut, RealType* restrict dist) const
  {
    return 0;
  }

  /// find index and distance of each nearest neighbor particle
  virtual void nearest_neighbor(std::vector<ripair>& ri,bool transposed=false) const
  {
    APP_ABORT("DistanceTableData::nearest_neighbor is not implemented in calling base class");
  }

  /// find indices and distances of nearest neighbors particles to particle n
  virtual void nearest_neighbors(int n,int neighbors,std::vector<ripair>& ri,bool transposed=false)
  {
    APP_ABORT("DistanceTableData::nearest_neighbors is not implemented in calling base class");
  }

  /// find species resolved indices and distances of nearest particles to particle n
  virtual void nearest_neighbors_by_spec(int n,int neighbors,int spec_start,std::vector<ripair>& ri,bool transposed=false)
  {
    APP_ABORT("DistanceTableData::nearest_neighbors is not implemented in calling base class");
  }

  inline void check_neighbor_size(std::vector<ripair>& ri,bool transposed=false) const
  {
    int m;
    if(transposed)
      m = N[SourceIndex];
    else
      m = N[VisitorIndex];
    if(ri.size()!=m)
      APP_ABORT("DistanceTableData::check_neighbor_size  distance/index vector length is not equal to the number of neighbor particles");
  }

  inline void print(std::ostream& os)
  {
    os << "Table " << Origin->getName() << std::endl;
    for(int i=0; i<r_m.size(); i++)
      os << r_m[i] << " ";
    os << std::endl;
  }

  const ParticleSet* Origin;

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
  /** full distance AB or AA  return r2_m(iat,jat) */
  Matrix<RealType,aligned_allocator<RealType> > r_m2;
  /** full displacement  AB or AA  */
  Matrix<PosType> dr_m2;
  /** J2 for compact neighbors */
  Matrix<int,aligned_allocator<int> > J2;
  /*@}*/

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
  void resize(int npairs, int nw)
  {
    N[WalkerIndex] =nw;
    //if(nw==1)
    {
      dr_m.resize(npairs);
      r_m.resize(npairs);
      rinv_m.resize(npairs);
      Temp.resize(N[SourceIndex]);
    }
  }

};
}
#endif
