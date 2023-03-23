//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: D. Das, University of Illinois at Urbana-Champaign
//                    Ken Esler, kpesler@gmail.com, University of Illinois at Urbana-Champaign
//                    Miguel Morales, moralessilva2@llnl.gov, Lawrence Livermore National Laboratory
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Raymond Clay III, j.k.rofling@gmail.com, Lawrence Livermore National Laboratory
//                    Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//                    Jeongnim Kim, jeongnim.kim@intel.com, Intel Corp
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////


#ifndef QMCPLUSPLUS_WALKER_H
#define QMCPLUSPLUS_WALKER_H

#include "OhmmsPETE/OhmmsMatrix.h"
#include "MinimalContainers/ConstantSizeMatrix.hpp"
#include "Pools/PooledData.h"
#include "Pools/PooledMemory.h"
#include "QMCDrivers/WalkerProperties.h"

namespace qmcplusplus
{
/** A container class to represent a walker.
 *
 * A walker stores the particle configurations {R}  and a property container.
 * RealTypehe template (P)articleSet(A)ttribute is a generic container  of position types.
 * RealTypehe template (G)radient(A)ttribute is a generic container of gradients types.
 * Data members for each walker
 * - ID : identity for a walker. default is 0.
 * - Age : generation after a move is accepted.
 * - Weight : weight to take the ensemble averages
 * - Multiplicity : multiplicity for branching. Probably can be removed.
 * - Properties  : 2D container. RealType first index corresponds to the H/Psi index and second index >=WP::NUMPROPERTIES.
 * - DataSet : a contiguous buffer providing a state snapshot of most/all walker data. 
     Much complicated state management arises in keeping this up to date, 
     or purposefully out of sync with actual datamembers. Here and in TWF, HAMs and PARTICLE sets
     associated with the walker.
 */
template<typename t_traits, typename p_traits>
class Walker
{
public:
  using WP = WalkerProperties::Indexes;
  enum
  {
    DIM = t_traits::DIM
  };
  /** typedef for real data type */
  using RealType = typename t_traits::RealType;
  /** typedef for estimator real data type */
  using FullPrecRealType = typename t_traits::FullPrecRealType;
  /** typedef for value data type. */
  using ValueType = typename t_traits::ValueType;
  /** array of particles */
  using ParticlePos = typename p_traits::ParticlePos;
  /** array of scalars */
  using ParticleScalar = typename p_traits::ParticleScalar;
  /** array of gradients */
  using ParticleGradient = typename p_traits::ParticleGradient;
  /** array of laplacians */
  using ParticleLaplacian = typename p_traits::ParticleLaplacian;
  /** typedef for value data type. */
  using SingleParticleValue = typename p_traits::SingleParticleValue;

  ///typedef for the property container, fixed size
  using PropertyContainer_t = ConstantSizeMatrix<FullPrecRealType, std::allocator<FullPrecRealType>>;

  /** @{
   * Not really "buffers", "walker message" also used to serialize walker, rename
   */
  using WFBuffer_t = PooledMemory<FullPrecRealType>;
  using Buffer_t   = PooledData<RealType>;
  /** }@ */

  ///id reserved for forward walking
  long ID;
  ///id reserved for forward walking
  long ParentID;
  ///DMCgeneration
  int Generation;
  ///Age of this walker age is incremented when a walker is not moved after a sweep
  int Age;
  ///Weight of the walker
  FullPrecRealType Weight;
  /** Number of copies for branching
   *
   * When Multiplicity = 0, this walker will be destroyed.
   */
  FullPrecRealType Multiplicity;
  /// mark true if this walker is being sent.
  bool SendInProgress;
  /// if true, this walker is either copied or tranferred from another MPI rank.
  bool wasTouched = true;

  /** The configuration vector (3N-dimensional vector to store
     the positions of all the particles for a single walker)*/
  ParticlePos R;

  //Dynamical spin variable.
  ParticleScalar spins;
#if !defined(SOA_MEMORY_OPTIMIZED)
  /** \f$ \nabla_i d\log \Psi for the i-th particle */
  ParticleGradient G;
  /** \f$ \nabla^2_i d\log \Psi for the i-th particle */
  ParticleLaplacian L;
#endif
  ///scalar properties of a walker
  PropertyContainer_t Properties;

  /** Property history vector
   *
   *  these are used as fixed length cyclic traces of a "property"
   */
  std::vector<std::vector<FullPrecRealType>> PropertyHistory;
  std::vector<int> PHindex;

  ///buffer for the data for particle-by-particle update
  WFBuffer_t DataSet;
  size_t block_end, scalar_end;

  // This is very useful for debugging transfer damage to walkers
#ifndef NDEBUG
private:
  bool has_been_on_wire_ = false;

public:
  bool get_has_been_on_wire() const { return has_been_on_wire_; }
  void set_has_been_on_wire(bool tf) { has_been_on_wire_ = tf; }
#endif

  ///create a walker for n-particles
  inline explicit Walker(int nptcl = 0)
      : Properties(1, WP::NUMPROPERTIES, 1, WP::MAXPROPERTIES)
  {
    ID                 = 0;
    ParentID           = 0;
    Generation         = 0;
    Age                = 0;
    Weight             = 1.0;
    Multiplicity       = 1.0;

    if (nptcl > 0)
      resize(nptcl);
    //static_cast<Matrix<FullPrecRealType>>(Properties) = 0.0;
  }

  Walker(const Walker& a) : Properties(1, WP::NUMPROPERTIES, 1, WP::MAXPROPERTIES) { makeCopy(a); }

  inline int addPropertyHistory(int leng)
  {
    int newL                            = PropertyHistory.size();
    PropertyHistory.push_back(std::vector<RealType>(leng, 0.0));
    PHindex.push_back(0);
    return newL;
  }

  inline void deletePropertyHistory() { PropertyHistory.clear(); }

  inline void resetPropertyHistory()
  {
    for (int i = 0; i < PropertyHistory.size(); i++)
    {
      PHindex[i] = 0;
      for (int k = 0; k < PropertyHistory[i].size(); k++)
      {
        PropertyHistory[i][k] = 0.0;
      }
    }
  }

  inline void addPropertyHistoryPoint(int index, FullPrecRealType data)
  {
    PropertyHistory[index][PHindex[index]] = (data);
    PHindex[index]++;
    if (PHindex[index] == PropertyHistory[index].size())
      PHindex[index] = 0;
    //       PropertyHistory[index].pop_back();
  }

  inline FullPrecRealType getPropertyHistorySum(int index, int endN)
  {
    FullPrecRealType mean = 0.0;
    typename std::vector<FullPrecRealType>::const_iterator phStart;
    phStart = PropertyHistory[index].begin() + PHindex[index];
    for (int i = 0; i < endN; phStart++, i++)
    {
      if (phStart >= PropertyHistory[index].end())
        phStart -= PropertyHistory[index].size();
      mean += (*phStart);
    }
    return mean;
  }

  ///assignment operator
  inline Walker& operator=(const Walker& a)
  {
    if (this != &a)
      makeCopy(a);
    return *this;
  }

  ///return the number of particles per walker
  inline int size() const { return R.size(); }

  ///resize for n particles
  inline void resize(int nptcl)
  {
    R.resize(nptcl);
    spins.resize(nptcl);
    G.resize(nptcl);
    L.resize(nptcl);
  }

  ///copy the content of a walker
  inline void makeCopy(const Walker& a)
  {
    ID                 = a.ID;
    ParentID           = a.ParentID;
    Generation         = a.Generation;
    Age                = a.Age;
    Weight             = a.Weight;
    Multiplicity       = a.Multiplicity;
    if (R.size() != a.R.size())
      resize(a.R.size());
    R = a.R;
    if (spins.size() != a.spins.size())
      resize(a.spins.size());
    spins = a.spins;
#if !defined(SOA_MEMORY_OPTIMIZED)
    G = a.G;
    L = a.L;
#endif
    Properties.copy(a.Properties);
    DataSet    = a.DataSet;
    block_end  = a.block_end;
    scalar_end = a.scalar_end;
    if (PropertyHistory.size() != a.PropertyHistory.size())
      PropertyHistory.resize(a.PropertyHistory.size());
    for (int i = 0; i < PropertyHistory.size(); i++)
      PropertyHistory[i] = a.PropertyHistory[i];
    PHindex = a.PHindex;
  }

  //return the address of the values of Hamiltonian terms
  inline FullPrecRealType* getPropertyBase() { return Properties.data(); }

  //return the address of the values of Hamiltonian terms
  inline const FullPrecRealType* getPropertyBase() const { return Properties.data(); }

  ///return the address of the i-th properties
  inline FullPrecRealType* getPropertyBase(int i) { return Properties[i]; }

  ///return the address of the i-th properties
  inline const FullPrecRealType* getPropertyBase(int i) const { return Properties[i]; }

  /** reset the property of a walker
   *@param logpsi \f$\log |\Psi|\f$
   *@param sigN  sign of the trial wavefunction
   *@param ene the local energy
   *
   *Assign the values and reset the age
   * but leave the weight and multiplicity
   */
  inline void resetProperty(FullPrecRealType logpsi, FullPrecRealType sigN, FullPrecRealType ene)
  {
    Age = 0;
    //Weight=1.0;
    Properties(WP::LOGPSI)      = logpsi;
    Properties(WP::SIGN)        = sigN;
    Properties(WP::LOCALENERGY) = ene;
  }

  /** reset the property of a walker
   * @param logpsi \f$\log |\Psi|\f$
   * @param sigN  sign of the trial wavefunction
   * @param ene the local energy
   * @param r2a \f$r^2\f$ for the accepted moves
   * @param r2p \f$r^2\f$ for the proposed moves
   * @param vq \f$\bar{V}/V\f$ scaling to control node divergency in JCP 93
   *
   *Assign the values and reset the age
   * but leave the weight and multiplicity
   */
  inline void resetProperty(FullPrecRealType logpsi,
                            FullPrecRealType sigN,
                            FullPrecRealType ene,
                            FullPrecRealType r2a,
                            FullPrecRealType r2p,
                            FullPrecRealType vq)
  {
    Age                         = 0;
    Properties(WP::LOGPSI)      = logpsi;
    Properties(WP::SIGN)        = sigN;
    Properties(WP::LOCALENERGY) = ene;
    Properties(WP::R2ACCEPTED)  = r2a;
    Properties(WP::R2PROPOSED)  = r2p;
    Properties(WP::DRIFTSCALE)  = vq;
  }

  /** marked to die
       *
       * Multiplicity and weight are set to zero.
       */
  inline void willDie()
  {
    Multiplicity = 0;
    Weight       = 0.0;
  }

  /** reset the walker weight, multiplicity and age */
  inline void reset()
  {
    Age          = 0;
    Multiplicity = 1.0e0;
    Weight       = 1.0e0;
  }

  inline void resizeProperty(int n, int m) { Properties.resize(n, m); }


  /** byte size for a packed message
   *
   * ID, Age, Properties, R, Drift, DataSet is packed
   */
  inline size_t byteSize()
  {
    // TODO: fix this! this is a non intuitive side effect for a size call
    //       breaks a bunch of things that could be const
    if (!DataSet.size())
    {
      registerData();
      DataSet.allocate();
    }
    return DataSet.byteSize();
  }

  void registerData()
  {
    // walker data must be placed at the beginning
    assert(DataSet.size() == 0);
    // scalars
    DataSet.add(ID);
    DataSet.add(ParentID);
    DataSet.add(Generation);
    DataSet.add(Age);
    // vectors
    assert(R.size() != 0);
    DataSet.add(R.first_address(), R.last_address());
    assert(spins.size() != 0);
    DataSet.add(spins.first_address(), spins.last_address());
#if !defined(SOA_MEMORY_OPTIMIZED)
    assert(G.size() != 0);
    DataSet.add(G.first_address(), G.last_address());
    assert(L.size() != 0);
    DataSet.add(L.first_address(), L.last_address());
#endif
    //Don't add the nLocal but the actual allocated size.  We want to register once for the life of a
    //walker so we leave space for additional properties.
    DataSet.add(Properties.data(), Properties.data() + Properties.capacity());
    //DataSet.add(Properties.first_address(), Properties.last_address());

    // \todo likely to be broken if the Properties change above is needed.
    for (int iat = 0; iat < PropertyHistory.size(); iat++)
      DataSet.add(PropertyHistory[iat].data(), PropertyHistory[iat].data() + PropertyHistory[iat].size());
    DataSet.add(PHindex.data(), PHindex.data() + PHindex.size());
    block_end  = DataSet.current();
    scalar_end = DataSet.current_scalar();
  }

  void copyFromBuffer()
  {
    assert(DataSet.size() != 0);
    DataSet.rewind();
    DataSet >> ID >> ParentID >> Generation >> Age;
    // vectors
    assert(R.size() != 0);
    DataSet.get(R.first_address(), R.last_address());
    assert(spins.size() != 0);
    DataSet.get(spins.first_address(), spins.last_address());
#if !defined(SOA_MEMORY_OPTIMIZED)
    assert(G.size() != 0);
    DataSet.get(G.first_address(), G.last_address());
    assert(L.size() != 0);
    DataSet.get(L.first_address(), L.last_address());
#endif
    DataSet.get(Properties.data(), Properties.data() + Properties.capacity());
    for (int iat = 0; iat < PropertyHistory.size(); iat++)
      DataSet.get(PropertyHistory[iat].data(), PropertyHistory[iat].data() + PropertyHistory[iat].size());
    DataSet.get(PHindex.data(), PHindex.data() + PHindex.size());
    assert(block_end == DataSet.current());
    assert(scalar_end == DataSet.current_scalar());
  }

  // In Clang 4.0.1 it is unclear why but attibute list cannot go at the end
  // when a template declaration follows
  void updateBuffer()
  {
    DataSet.rewind();
    DataSet << ID << ParentID << Generation << Age;
    // vectors
    DataSet.put(R.first_address(), R.last_address());
    DataSet.put(spins.first_address(), spins.last_address());
#if !defined(SOA_MEMORY_OPTIMIZED)
    DataSet.put(G.first_address(), G.last_address());
    DataSet.put(L.first_address(), L.last_address());
#endif
    DataSet.put(Properties.data(), Properties.data() + Properties.capacity());
    for (int iat = 0; iat < PropertyHistory.size(); iat++)
      DataSet.put(PropertyHistory[iat].data(), PropertyHistory[iat].data() + PropertyHistory[iat].size());
    DataSet.put(PHindex.data(), PHindex.data() + PHindex.size());
    assert(block_end == DataSet.current());
    assert(scalar_end == DataSet.current_scalar());
  }

  template<class Msg>
  inline Msg& putMessage(Msg& m)
  {
    // Pack DataSet buffer
    if (!DataSet.size())
    {
      registerData();
      DataSet.allocate();
    }
    updateBuffer();
    m.Pack(DataSet.data(), DataSet.size());
    return m;
  }

  template<class Msg>
  inline Msg& getMessage(Msg& m)
  {
    if (!DataSet.size())
    {
      registerData();
      DataSet.allocate();
    }
    m.Unpack(DataSet.data(), DataSet.size());
    copyFromBuffer();
    return m;
  }
};

template<class RealType, class PA>
std::ostream& operator<<(std::ostream& out, const Walker<RealType, PA>& rhs)
{
  copy(rhs.Properties.begin(), rhs.Properties.end(), std::ostream_iterator<double>(out, " "));
  out << std::endl;
  out << rhs.R;
  return out;
}
} // namespace qmcplusplus

#endif
