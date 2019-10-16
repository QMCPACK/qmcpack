//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by:
//
// File created by: Jeongnim Kim, jeongnim.kim@intel.com, Intel Corp.
//////////////////////////////////////////////////////////////////////////////////////


/**@file ParticleAttrib.h
 *
 * Declaraton of ParticleAttrib<T> derived from Vector
 */

#ifndef OHMMS_NEW_PARTICLEATTRIB_PEPE_H
#define OHMMS_NEW_PARTICLEATTRIB_PEPE_H

#include <OhmmsPETE/OhmmsVector.h>
#include "Utilities/OhmmsObject.h"
#include "Lattice/CrystalLattice.h"

namespace qmcplusplus
{

/** Attaches a unit to a Vector for IO
 *
 *  Makes Vect
 */
template<class T, typename Alloc = std::allocator<T>>
class ParticleAttrib : public Vector<T, Alloc>, public OhmmsObject
{
  typedef Vector<T, Alloc> __my_base;

public:
  /// The unit type
  PosUnit InUnit;

  /** constructor with size n*/
  explicit inline ParticleAttrib(size_t n = 0) : __my_base(n), InUnit(PosUnit::Cartesian) {}

  /** constructor with an initialized ref */
  explicit inline ParticleAttrib(T* ref, size_t n) : __my_base(ref, n), InUnit(PosUnit::Cartesian) {}

  ParticleAttrib(const ParticleAttrib& rhs) = default;
  inline ParticleAttrib& operator=(const ParticleAttrib& rhs) = default;

  /** assignment operator to enable PETE */
  template<class RHS>
  inline ParticleAttrib& operator=(const RHS& rhs)
  {
    assign(*this, rhs);
    return *this;
  }

  //@{set/set the unit
  inline void setUnit(PosUnit i) { InUnit = i; }
  inline PosUnit getUnit() const { return InUnit; }
  //@}

  OhmmsObject* makeClone() const { return new ParticleAttrib<T, Alloc>(*this); }

  /** Specialized to write the unit
   *\return true, if the attribute is relative to a unit
   *
   * Ad hoc get function to tell users if this attribute is absolute or relative value.
   * E.g., is Cartesian unit vs Lattice unit.
   */
  bool get(std::ostream& os) const
  {
    os << InUnit;
    return true;
  }

  /** Specialized to read the unit
   */
  bool put(std::istream& is)
  {
    is >> InUnit;
    return true;
  }

  /*@warning not fully implemented.*/
  bool put(xmlNodePtr cur) { return true; }

  ///reset member data
  void reset() {}
};

template<class T, unsigned D>
inline T* get_first_address(ParticleAttrib<TinyVector<T, D>>& a)
{
  return &(a[0][0]);
}

template<class T, unsigned D>
inline T* get_last_address(ParticleAttrib<TinyVector<T, D>>& a)
{
  return &(a[0][0]) + D * a.size();
}

} // namespace qmcplusplus

#include "ParticleBase/ParticleAttrib.cpp"

#endif // OHMMS_PARTICLEATTRIB_PEPE_H
