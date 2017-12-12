//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////
    
    



#ifndef QMCPLUSPLUS_RADIALGRIDFUNCTOR_BUILDERBASE_H
#define QMCPLUSPLUS_RADIALGRIDFUNCTOR_BUILDERBASE_H

#include "Configuration.h"
#include "Numerics/OneDimGridFunctor.h"
#include "QMCWaveFunctions/SphericalOrbitalSet.h"
#include "OhmmsData/HDFAttribIO.h"

namespace qmcplusplus
{

/**An abstract base class to manage the radial grids and orbitals that are associated with a center.
 *
 */
struct RGFBuilderBase: public QMCTraits
{

  typedef OneDimGridBase<RealType>                        GridType;
  typedef OneDimGridFunctor<RealType>                     RadialOrbitalType;
  typedef SphericalOrbitalSet<RadialOrbitalType,GridType> CenteredOrbitalType;
  bool debug;

  ///the species
  std::string m_species;
  ///the radial orbitals
  CenteredOrbitalType* m_orbitals;

  RGFBuilderBase(): m_orbitals(NULL), Counter(0), m_file_id(-1), m_group_id(-1),debug(false) { }

  virtual ~RGFBuilderBase();

  ///assign a CenteredOrbitalType to work on
  void setOrbitalSet(CenteredOrbitalType* oset, const std::string& acenter)
  {
    m_orbitals = oset;
    m_species = acenter;
  }

  ///add a radial grid
  virtual bool addGrid(xmlNodePtr cur);

  ///add a radial orbital with quantum numbers nlms
  virtual bool addRadialOrbital(xmlNodePtr cur, const QuantumNumberType& nlms) = 0;

  ///add common parameters
  virtual bool putCommon(xmlNodePtr cur) = 0;

  void print(const std::string& aroot, int omode, bool debug);

  ///keeps track of number of times HDF5 file has been accessed
  int Counter;
  ///id for HDF5 file
  hid_t m_file_id;
  ///id for main group of HDF5 file
  hid_t m_group_id;
};
}
#endif
