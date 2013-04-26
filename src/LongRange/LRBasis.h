#ifndef QMCPLUSPLUS_LRBASIS_H
#define QMCPLUSPLUS_LRBASIS_H

#include "Particle/ParticleSet.h"
#include "Utilities/PooledData.h"

namespace qmcplusplus
{

/** @ingroup longrange
 *\brief Base-class for long-range breakups.
 *
 * Contains 3 important functions: c(n,k), h(n,r), hintr2(r)
 *  which evaluate the n'th basis function at k or r.
 *  \f[ hintr2_n = \int dr h_n(r) r^2, \f]
 */

class LRBasis: public QMCTraits
{
protected:
  ///size of the basis elements
  int BasisSize;
  ///Real-space cutoff for short-range part
  RealType m_rc;

  ///Typedef for the lattice-type. We don't need the full particle-set.
  typedef ParticleSet::ParticleLayout_t ParticleLayout_t;

  //A copy of the lattice, so that we have the cell-vectors + volume.
  ParticleLayout_t& Lattice;

public:

  //Constructor: copy lattice and init rc
  LRBasis(ParticleLayout_t& ref): Lattice(ref), m_rc(0.0)
  {
    /*Do nothing*/
  }

  virtual ~LRBasis() { }

  inline int NumBasisElem() const
  {
    return BasisSize;
  }

  //Real-space basis function + integral: override these
  virtual RealType h(int n, RealType r) const = 0;
  virtual RealType hintr2(int n) = 0;
  //k-space basis function: override this
  virtual RealType c(int m, RealType k) = 0;

  //May need extra functionality when resetting rc. Override this.
  virtual void set_rc(RealType rc) = 0;
  inline RealType get_rc()
  {
    return m_rc;
  }
  inline RealType get_CellVolume()
  {
    return Lattice.Volume;
  }
  inline ParticleLayout_t& get_Lattice()
  {
    return Lattice;
  }
  inline void set_Lattice(ParticleLayout_t& ref)
  {
    Lattice = ref;
  }
};

}

#endif
